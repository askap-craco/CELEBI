#!/usr/bin/env python
"""
Abstracts away VCRAFT files

Copyright (C) CSIRO 2017
"""
__author__ = "Keith Bannister <keith.bannister@csiro.au>"

import itertools
import logging
import os
import sys
import warnings

import freqconfig
import numpy as np
from crafthdr import DadaHeader

log = logging.getLogger(__name__)

bat_cards = ("START_WRITE_BAT40", "STOP_WRITE_BAT40", "TRIGGER_BAT40")
frame_cards = ("START_WRITE_FRAMEID", "STOP_WRITE_FRAMEID", "TRIGGER_FRAMEID")

# number of samples per 32 bit word indexed by mode
# see setupCraftDownload
# https://svn.atnf.csiro.au/askapsoft/Src/trunk/Code/Base/adbe/current/adbe/Redback3.cc
#
# Andrew's comment from setupCraftDownload below:
# // For a DIMM with 28 address bits there are (2^(28-7))/72=29172 sets of 72 beams
# // The 28-7 comes from the fact that the lower 2 bits are used for address incrementing
# // and the next 5 bits are for counting samples (up to 32). So the remaining number
# // of bits left for sets of 72 beams is 28-7 = 21 bits. We divide by 72 because
# // that's how many beams we have in a set. Each set contains 32 samples for each beam
# // (assuming 16bit data). Therefore the max number of timeSamples that can be read is
# // 29172 sets x 32 samples. For single beam mode we have 36 more samples available but
# // only for the 2 beams selected

SAMPS_PER_WORD32 = [1, 2, 4, 16, 1, 2, 4, 16]
MODE_BEAMS = [72, 72, 72, 72, 2, 2, 2, 2]
MAX_SETS = 29172

#  Number of sample vs mode
SAMPS_MODE = [
    29172 * 32 * nsamp_per_word * 72 / nbeams
    for (nsamp_per_word, nbeams) in zip(SAMPS_PER_WORD32, MODE_BEAMS)
]


def unpack_craft(fin, nsamp_per_word):
    dwords = np.fromfile(fin, dtype=np.uint32)
    nwords = len(dwords) / nchan
    nsamps = nwords * nsamp_per_word
    assert len(dwords) == nchan * nwords
    dwords.shape = (nwords, nchan)

    d = np.empty((nsamps, nchan, 2), dtype=np.int8)
    bitshift = 32 / nsamp_per_word
    assert nsamp_per_word * bitshift == 32

    for samp in range(nsamp_per_word):
        # convert from 4-bit two's complement to int8
        d[samp::4, :, 0] = (dwords & 0xF) - (dwords & 0x8) * 2  # real
        dwords >> 4
        d[samp::4, :, 1] = (dwords & 0xF) - (dwords & 0x8) * 2  # imag
        dwords >> 4


class VcraftFile:
    def __init__(self, fname, mode=None):
        self.fname = fname
        f = fname
        hdr = DadaHeader.fromfile(
            f + ".hdr"
        )  # hdr is a dictionary of the header contents
        self.hdr = hdr
        self.fin = open(f)
        self.hdrsize = int(hdr["HDR_SIZE"][0])
        self.fin.seek(
            self.hdrsize
        )  # Set the file's position to after the header
        if mode is None:
            self.mode = int(hdr["CRAFT_MODE"][0]) & 0x3
        else:
            self.mode = mode

        self.nbits = int(hdr["NBITS"][0])
        self.infsamp = float(hdr["SAMP_RATE"][0])  # samples per second
        # TODO: Calculate frequencies a bit better - tricky because individual
        # files have freqs with gaps, which makes life a little wierd in sigprocland
        self.freqs = np.array(list(map(float, hdr["FREQS"][0].split(","))))

        # TODO: Get tstart from BATs. This is the easy way
        # self.tstart = float(hdr['ANT_MJD'][0])
        # Trigger is at the *end* of the buffer.
        self.trigger_frameid = int(hdr["TRIGGER_FRAMEID"][0])
        self.trigger_mjd = float(hdr["TRIGGER_MJD"][0])
        self.expected_nsamp = int(
            hdr.get_value("NSAMPS_REQUEST", SAMPS_MODE[self.mode])
        )
        self.start_frameid = self.trigger_frameid - self.expected_nsamp
        self.start_mjd = (
            self.trigger_mjd - self.expected_nsamp / self.infsamp / 86400.0
        )
        self.bats = np.array([int(hdr[c][0], base=16) for c in bat_cards])
        self.frames = np.array([int(hdr[c][0]) for c in frame_cards])
        self.bat0 = int(hdr["NOW_BAT"][0], base=16)
        self.bat0_40 = self.bat0 & 0xFFFFFFFFFF
        self.beam = int(hdr["BEAM"][0])
        self.pol = hdr.get("POL", (None, None))[0]
        if self.pol is None:
            if self.beam % 2 == 0:
                self.pol = "X"
            else:
                self.pol = "Y"

        if self.nsamps != self.expected_nsamp:
            warnings.warn(
                f"VCRAFT file {self.fname} had unexpected size. Mode {self.mode} expected {self.expected_nsamp} actual {self.nsamps}"
            )

    def __str__(self):
        return f"{self.fname} NBITS={self.nbits}"

    __repr__ = __str__

    @property
    def nsamps(self):
        """Return the number of complex samples in the file"""

        file_bytes = os.path.getsize(self.fname)
        data_bytes = file_bytes - self.hdrsize
        nchan = len(self.freqs)
        smp_per_word = SAMPS_PER_WORD32[self.mode]
        data_words = data_bytes / 4

        nsamp = data_words * smp_per_word / nchan

        return nsamp

    def print_times(self):
        bats = self.bats
        frames = self.frames

        print(
            "BAT duration (s)",
            (bats[0] - bats[1]) / 1e6,
            "offset",
            (bats[2] - bats[0]) / 1e6,
            "file offset",
            (bat0_40 - bats[0]) / 1e6,
        )
        # lowest 32 bits of bat from the time the file was written
        print(
            "FRAMES duration",
            (frames[0] - frames[1]) / infsamp,
            "offset",
            (frames[2] - frames[0]) / infsamp,
        )

    def read(self, startsamp=0, nsamp=None):
        """Reads everything and returns a numpy array
        with shape dtype=(nsamps,nchan) and dtype=np.complex64
        """
        # See https://svn.atnf.csiro.au/askapsoft/Src/trunk/Code/Base/adbe/current/adbe/Beamformer.cc
        mode = self.mode
        fin = self.fin
        nchan = len(self.freqs)
        self.fin.seek(self.hdrsize)
        assert startsamp >= 0, f"Invalid startsamp {startsamp} {self.fname}"
        if nsamp is None:
            nsamp = self.nsamps - startsamp

        assert (
            startsamp + nsamp <= self.nsamps
        ), f"Asked for too many samples. startsamp {startsamp} nsamp {nsamp} nsamps {self.nsamps}"

        if mode == 0:  # 16b+16b
            offset_bytes = startsamp * nchan * 4  # each sample is 4 bytes
            self.fin.seek(self.hdrsize + offset_bytes)
            d = np.fromfile(fin, dtype=np.int16, count=nsamp * nchan * 2)
            # truncate if required
            # d = d[:2*nchan*nsamps]
            assert 2 * nchan * nsamp == len(
                d
            ), "Not integral number of samples"
            d.shape = (nsamp, nchan, 2)
            nsamps = nsamp
        elif mode == 1:  # 8b + 8b:
            # This works perfectly, but we need some practice at the other method with non-native types
            # d = np.fromfile(fin, dtype=np.int8, count=nsamp*nchan*2) # 2 for complex
            # nsamps = len(d)/nchan/2
            # assert 2*nchan*nsamps == len(d), 'Not integral number of samples'
            # print 'MODE1', startsamp, 'requested ', nsamp, 'got', nsamps, 'count', nsamp*nchan*2, 'len', len(d)
            # each 32 bit word contains 4x8 bit numbers imag/real/imag/real is two seuential samplesfrom a single channel
            # d.shape = (nsamps, nchan, 2, 2)
            # let's put those axes next to each other
            # d = d.transpose(0,2,1,3)
            # and reshape it to the orrect shape
            # d = d.reshape(nsamps, nchan, 2)

            wordidx = startsamp / 2  # which 32 bit word the start sample is in
            sampoff = (
                startsamp % 2
            )  # how may samples into the first word the start sample is
            nwordsamps = (
                nsamp + 1 + sampoff
            ) / 2  # how many times we need to read (in words)
            nwords = (
                nwordsamps * nchan
            )  # total numberof words including channels
            seek_bytes = (
                self.hdrsize + wordidx * nchan * 4
            )  # seek offset in bytes
            fin.seek(seek_bytes)
            # print 'MODE1', 'startsamp', startsamp, 'wordidx', wordidx, 'sampoff', sampoff, 'nwordsamps', nwordsamps, 'nwords', nwords, 'seek bytes', seek_bytes
            dwords = np.fromfile(fin, dtype="<u4", count=nwords)

            # each word contains 4 8 bit numbers, (imag/real)*2
            nwords = len(dwords) / nchan
            assert (
                len(dwords) == nchan * nwords
            ), f"Got {len(dwords)} dwords = nchan={nchan} nwords={nwords} expected={nchan*nwords}"
            dwords.shape = nwords, nchan
            nsamps = nwords * 2
            d = np.empty((nsamps, nchan, 2), dtype=np.int8)
            for samp in range(2):
                # convert from 8-bit two's complement to int8
                d[samp::2, :, 0] = (dwords & 0xFF) - (
                    dwords & 0x80
                ) * 2  # real
                dwords >>= 8
                d[samp::2, :, 1] = (dwords & 0xFF) - (
                    dwords & 0x80
                ) * 2  # imag
                dwords >>= 8

            d = d[sampoff : sampoff + nsamp, :, :]
            nsamps = nsamp

        # AT LEAST FRB200430 ak01 uses mode 2 (do all CRAFT FRBs?)
        elif mode == 2:  # 4b+4b

            wordidx = startsamp / 4  # which 32 bit word the start sample is in
            sampoff = (
                startsamp % 4
            )  # how may samples into the first word the start sample is
            nwordsamps = (
                nsamp + 3 + sampoff
            ) / 4  # how many times we need to read (in words)
            nwords = (
                nwordsamps * nchan
            )  # total numberof words including channels
            seek_bytes = (
                self.hdrsize + wordidx * nchan * 4
            )  # seek offset in bytes
            fin.seek(seek_bytes)
            dwords = np.fromfile(fin, dtype="<u4", count=nwords)
            # each word contains 4 8 bit numbers, (imag/real)*2
            nwords = len(dwords) / nchan
            nsamps = nwords * 4
            assert len(dwords) == nchan * nwords
            dwords.shape = (nwords, nchan)
            d = np.empty((nsamps, nchan, 2), dtype=np.int8)
            nsamps = nsamp
            for samp in range(4):
                # convert from 4-bit two's complement to int8
                d[samp::4, :, 0] = (dwords & 0xF) - (dwords & 0x8) * 2  # real
                dwords >>= 4
                d[samp::4, :, 1] = (dwords & 0xF) - (dwords & 0x8) * 2  # imag
                dwords >>= 4

            d = d[sampoff : sampoff + nsamp, :, :]

        elif mode == 3:  # 1b+1b
            # each 32 bit word contais 32 1 bit numbers (imag/real)*16 for the same channel
            wordidx = (
                startsamp / 16
            )  # which 32 bit word the start sample is in
            sampoff = (
                startsamp % 16
            )  # how may samples into the first word the start sample is
            nwordsamps = (
                nsamp + 15 + sampoff
            ) / 16  # how many times we need to read (in words)
            nwords = (
                nwordsamps * nchan
            )  # total numberof words including channels
            seek_bytes = (
                self.hdrsize + wordidx * nchan * 4
            )  # seek offset in bytes
            fin.seek(seek_bytes)
            dwords = np.fromfile(fin, dtype="<u4", count=nwords)
            # nwords = len(dwords)/nchan
            assert len(dwords) == nwords
            dwords.shape = (nwordsamps, nchan)
            nsamps = nsamp
            d = np.empty((nwordsamps * 16, nchan, 2), dtype=np.int8)
            for samp in range(16):
                # Convert 0 into +1 and 1 into -1
                d[samp::16, :, 0] = 1 - 2 * (dwords & 0x1)  # real
                dwords >>= 1
                d[samp::16, :, 1] = 1 - 2 * (dwords & 0x1)  # imag
                dwords >>= 1
            log.debug(
                "MODE3 startsamp=%s sampoff=%s wordidx=%s nwordssamps=%s nwords=%s seek_bytes=%s nsamps=%s dwords.shape=%s d.shape=%s",
                startsamp,
                sampoff,
                wordidx,
                nwordsamps,
                nwords,
                seek_bytes,
                nsamps,
                dwords.shape,
                d.shape,
            )
            d = d[sampoff : sampoff + nsamp, :, :]
            assert (
                d.shape[0] == nsamp
            ), f"Incorrect output shape {d.shape} expected {nsamp}"
        else:
            raise ValueError(f"Unsupported mode (yet) {mode}")

        df = np.empty((nsamps, nchan), dtype=np.complex64)
        df.real = d[:, :, 0]
        df.imag = d[:, :, 1]
        return df

    def print_summary(self):
        bats = np.array([int(self.hdr[b][0], base=16) for b in bat_cards])
        frames = np.array([int(self.hdr[b][0], base=10) for b in frame_cards])
        for b, bc in zip(bat_cards, bats):
            print(
                b, self.hdr[b][0], bc, bc - bats[0], (bc - bats[0]) / 1e6, "s"
            )

        for b, bc in zip(frame_cards, frames):
            print(
                b,
                self.hdr[b][0],
                bc,
                bc - frames[0],
                (bc - frames[0]) / (1e6 * 32 / 27),
                "s",
            )


class VcraftMux:
    """
    Multiplexes together VCRAFT files to give a contiguous, monotonically increasing
    frequency axis
    """

    # Multiplex: a system or signal involving simultaneous transmission of several messages along a single channel of communication

    def __init__(self, vcraft_files, delays=None):
        """
        :vcraft_files: A list of open Vcraft files
        """
        self._files = sorted(
            vcraft_files,
            key=lambda f: (f.hdr["CARD_NO"][0], f.hdr["FPGA_ID"][0]),
        )
        if delays is None:
            delays = {}

        self.file_delays = np.array(
            [
                int(delays.get(os.path.basename(f.fname + ".hdr"), 0))
                for f in vcraft_files
            ]
        )
        self.ant = self.hdr_identical("ANT")
        self.antno = int(self.hdr_identical("ANTENNA_NO"))
        self.beam = int(self.hdr_identical("BEAM"))
        self.nbits = int(self.hdr_identical("NBITS"))
        self.mode = int(self.hdr_identical("MODE"))
        self.data_type = self.hdr_identical("DATA_TYPE")  # just a check
        assert self.data_type == "CRAFT_VOLTAGES"
        self.samp_rate = float(self.hdr_identical("SAMP_RATE"))
        freq_str = self.allhdr("FREQS")
        freqs = np.array(
            [list(map(float, flist.split(","))) for flist in freq_str]
        )

        # there was traditionally a frequency offset of + 1inserted in the header
        # See craft-232
        # undo it if nothing is specified in the header to the contrary

        freq_offset = float(self.hdr_identical("FREQ_OFFSET", -1.0))
        freqs += freq_offset

        nchan_per_file = len(freqs[0])
        assert freqs.shape == (len(self._files), nchan_per_file)

        self.freqconfig = freqconfig.FreqConfig(freqs, reverse=True)
        self.freqs = (
            np.arange(self.freqconfig.nchan_span) * self.freqconfig.bw
            + self.freqconfig.freq
        )

        self.all_samps = [f.nsamps for f in self._files]
        self.nsamps = min(self.all_samps)

        self.trigger_frameids = np.array(
            list(map(int, self.allhdr("TRIGGER_FRAMEID")))
        )
        self.trigger_mjds = np.array(
            list(map(float, self.allhdr("TRIGGER_MJD")))
        )

        self.start_mjds = np.array([f.start_mjd for f in vcraft_files])
        self.start_mjd = max(self.start_mjds)

        self.start_frameids = (
            np.array([f.start_frameid for f in vcraft_files])
            - self.file_delays
        )
        self.start_frameid = max(self.start_frameids)
        # BUG! Start MJDs don't account for sample offsets
        self.sample_offsets = (
            self.start_frameid
            - np.array([f.start_frameid for f in self._files])
            + np.array(
                [
                    int(delays.get(os.path.basename(f.fname + ".hdr"), 0))
                    for f in self._files
                ]
            )
        )  # temporary
        # self.sample_offsets = self.start_frameid - self.start_frameids

        self.all_overlap_samps = [
            f.nsamps - self.sample_offsets[i]
            for i, f in enumerate(self._files)
        ]
        self.overlap_nsamps = min(self.all_overlap_samps)

        print(
            "SAMPLE OFFSETS",
            self.sample_offsets,
            "FILE DELAYS",
            self.file_delays,
        )
        assert np.all(self.sample_offsets >= 0)
        assert np.all(self.sample_offsets < self.nsamps)
        self.start_mjd = self.start_mjds[np.argmin(self.sample_offsets)]
        assert np.all(
            abs(self.start_mjds - self.trigger_mjds) < 1
        ), "MJD adjustment should be << 1 day"

        beam_ra = np.array(
            list(map(float, self.allhdr("BEAM_RA")))
        ).mean()  # TODO: make craft_vdump write same BEAM_RA for all card/fpgas
        beam_dec = np.array(list(map(float, self.allhdr("BEAM_DEC")))).mean()
        # self.beam_pos = SkyCoord(beam_ra, beam_dec, frame='icrs', unit=('deg','deg'))
        self.beam_pos = (beam_ra, beam_dec)
        self.hdr = self._files[0].hdr
        pols = np.array([f.pol for f in self._files])
        assert np.all(pols == pols[0])
        self.pol = pols[0]

    def allhdr(self, cardname, default=None):
        """
        Returns a list of header values for all files, with the given card name
        """
        header_values = [
            f.hdr.get_value(cardname, default) for f in self._files
        ]
        return header_values

    def hdr_identical(self, cardname, default=None):
        """
        Returns the header value for the given cardname.
        Checks the value is the same for all files. If not, it throws a ValueError
        """
        hdr_values = set(self.allhdr(cardname, default))
        if len(hdr_values) != 1:
            raise ValueError(
                f"Exepected the same header value for {cardname} for all files. Got these values {self.allhdr(cardname, default)}"
            )

        return hdr_values.pop()

    def read(self, samp_start=0, nsamp=None):
        """
        Read into a giant buffer - demuxing along the way.
        Probably not memory optimal, but it'll have to do for now.
        Something with 'yield' in it would probably make more sense - I should do that
        """
        if nsamp is None:
            nsamp = self.nsamps

        assert (
            samp_start + nsamp <= self.nsamps
        ), f"Invalid read request. User gave: samp_start={samp_start}, nsamp={nsamp}. File nsmaps={self.nsamps}"
        print(f"#!#!# n={nsamp} is VALID!")

        # allocate giant buffer
        d = np.zeros((nsamp, len(self.freqs)), dtype=np.complex64)
        for ifile, f in enumerate(self._files):
            out_chans = self.freqconfig.chanmaps[ifile, :]
            fsamp_start = samp_start + self.sample_offsets[ifile]
            try:
                d[:, out_chans] = f.read(fsamp_start, nsamp)
            except:
                warnings.warn("File read at funny time. Setting to zero")

        return d


def mux_by_pol(filenames, delays=None):
    """
    :return: Dictionary keyened by 'X' or "Y
    """
    all_files = [VcraftFile(f) for f in filenames]
    all_files.sort(key=lambda f: f.pol)
    mux_by_pol = itertools.groupby(all_files, lambda f: f.pol)
    muxes = {pol: VcraftMux(list(files), delays) for pol, files in mux_by_pol}

    return muxes


def mux_by_antenna(filenames, delays=None):
    all_files = [VcraftFile(f) for f in filenames]
    all_files.sort(key=lambda f: f.hdr["ANT"][0])

    # ants = [f.hdr['ANT'][0] for f in all_files]     # Redundant line

    mux_by_ant = itertools.groupby(
        all_files, lambda f: f.hdr["ANT"][0]
    )  # An iterable of tuples (antname, files) grouped by antenna name

    muxes = [VcraftMux(list(files), delays) for antname, files in mux_by_ant]
    muxes.sort(key=lambda mux: mux.antno)  # sort by antenna number

    return muxes


def _main():
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

    parser = ArgumentParser(
        description="Script description",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Be verbose",
    )
    parser.add_argument(
        "-n",
        "--nsamps",
        help="Number of samples per integration",
        type=int,
        default=1500,
    )
    parser.add_argument(
        "-s", "--show", help="Show plots", action="store_true", default=False
    )
    parser.add_argument(dest="files", nargs="+")
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    all_files = []
    for f in values.files:
        print(f)
        vf = VcraftFile(f)
        vf.print_summary()
        all_files.append(vf)

    muxes = mux_by_antenna(values.files)
    for mux in muxes:
        print(mux.freqconfig)
        print(mux.freqconfig.freqs)
        print(mux.freqconfig.chanmaps)
        print(mux.freqconfig.freqmaps)
        print(mux.freqs)


if __name__ == "__main__":
    _main()
