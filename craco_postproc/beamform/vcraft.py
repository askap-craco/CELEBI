import numpy as np
import itertools
import freqconfig

from pathlib import Path
from astropy.time import Time, TimeDelta
from craftpy.crafthdr import CraftHeader
from sigpyproc import libSigPyProc


class UnpackInfo(object):
    mode_to_nbits = {0: 16, 1: 8, 2: 4, 3: 1}
    # 16b+16b, 8b+8b, 4b+4b, 1b+1b

    nbits_to_dtype = {
        1: "<u1",
        2: "<u1",
        4: "<u1",
        8: "<u1",
        16: "<u2",
        32: "<f4",
    }

    def __init__(self, mode):
        self._mode = mode

    @property
    def mode(self):
        return self._mode

    @property
    def nbits(self):
        """Number of bits per component (real/Imag)"""
        return self.mode_to_nbits[self.mode & 0x3]

    @property
    def dtype(self) -> np.dtype:
        return np.dtype(self.nbits_to_dtype[self.nbits])

    @property
    def itemsize(self) -> int:
        return self.dtype.itemsize

    @property
    def unpack(self) -> bool:
        return bool(self.nbits in {1, 2, 4})

    @property
    def bitfact(self) -> int:
        return 8 // self.nbits if self.unpack else 1

    @property
    def sampsize(self):
        return self.itemsize / self.bitfact


class VCRAFTHeader(CraftHeader):
    """
    To read CRAFT antenna header files
    return: A ordered dict of all header flags including beam positions.
    """

    @property
    def hdr_size(self):
        return int(self.HDR_SIZE)

    @property
    def nbits(self):
        return int(self.NBITS)

    @property
    def beam(self):
        return int(self.BEAM)

    @property
    def ant(self):
        return str(self.ANT)

    @property
    def card_no(self):
        return int(self.CARD_NO)

    @property
    def fpga_id(self):
        return int(self.FPGA_ID)

    @property
    def pol(self):
        if self.POL is None:
            return "Y" if self.beam % 2 else "X"
        return self.POL

    @property
    def freqs(self):
        return np.array(self.FREQS.split(","), dtype=float)

    @property
    def nchans(self):
        return self.freqs.size

    @property
    def beam_ra(self):
        return float(self.BEAM_RA)

    @property
    def beam_dec(self):
        return float(self.BEAM_DEC)

    @property
    def craft_mode(self):
        return int(self.CRAFT_MODE)

    @property
    def mode(self):
        return int(self.MODE)

    @property
    def samp_rate(self):
        return float(self.SAMP_RATE)

    @property
    def nsamps_request(self):
        return int(self.NSAMPS_REQUEST) if self.NSAMPS_REQUEST else None

    @property
    def start_write_mjd(self):
        return float(self.START_WRITE_MJD)

    @property
    def stop_write_mjd(self):
        return float(self.STOP_WRITE_MJD)

    @property
    def trigger_mjd(self):
        return float(self.TRIGGER_MJD)

    @property
    def start_mjd(self):
        trigger = Time(self.trigger_mjd, format="mjd", scale="utc")
        start = trigger - TimeDelta(self.nsamps / self.samp_rate, format="sec")
        return start.mjd

    @property
    def trigger_bat(self):
        return str(self.TRIGGER_BAT)

    @property
    def trigger_frameid(self):
        return int(self.TRIGGER_FRAMEID)

    @property
    def start_frameid(self):
        return self.trigger_frameid - self.nsamps

    @property
    def start_bat(self):
        return int(self.trigger_bat - (self.trigger_frameid * 27 / 32))

    @property
    def bat0(self):
        return int(self.start_bat + 999999)

    @property
    def frame0(self):
        return int((self.bat0 - self.start_bat) * 32 / 27)


class VcraftFile(object):
    def __init__(self, filename):
        self._filename = filename
        self._header = VCRAFTHeader.fromfile(self.header_file)
        self._unpackinfo = UnpackInfo(self.mode)

    def __str__(self):
        return f"{self.filename} NBITS={self.header.nbits}"

    __repr__ = __str__

    @property
    def filename(self) -> str:
        """Name of the input file (`str`, read-only)."""
        return self._filename

    @property
    def header_file(self):
        """Header file (`str`, read-only)."""
        hdrfile = f"{self.filename}.hdr"
        if Path(hdrfile).is_file():
            return hdrfile
        else:
            raise IOError(f"Header file '{hdrfile}' does not exist.")

    @property
    def header(self):
        """Header metadata of input file (`str`, read-only)."""
        return self._header

    @property
    def mode(self):
        return self.header.craft_mode & 0x3

    @property
    def unpackinfo(self):
        return self._unpackinfo

    @property
    def nsamps(self):
        """ Return the number of complex samples in the file"""
        file_size = Path(self.filename).stat().st_size
        data_bytes = file_size - self.header.hdrsize
        return 8 * data_bytes // self.header.nbits // self.header.nchans

    def read(self, start=0, nsamps=None):
        """Reads everything and returns a numpy array
        with shape dtype=(nchan, nsamps) and dtype=np.complex64
        """
        if nsamps is None:
            nsamps = self.nsamps - start

        assert start >= 0, f"Invalid start samp to read: {start}"
        assert nsamps > 0, f"Invalid samples to read: {nsamps}, start: {start}"
        assert (
            start + nsamps <= self.nsamps
        ), f"Asked for too many samples. start samp: {start} nsamps: {nsamps} nsamps in file: {self.nsamps}"

        count = int(self.header.nchans * nsamps // self.bitfact) * 2
        data = np.fromfile(
            self.filename,
            count=count,
            dtype=self.unpackinfo.dtype,
            offset=self.header.hdrsize
            + start * self.unpackinfo.sampsize * 2 * self.header.nchans,
        )
        if self.unpackinfo.unpack:
            data = libSigPyProc.unpack(data, self.unpackinfo.nbits)

        nsamps_read = data.size // self.header.nchans // 2
        dfout = data.astype(np.float32).view(np.complex64)
        dfout = dfout.reshape(nsamps_read, self.header.nchans).transpose()
        return dfout


class VcraftMux(object):
    """
    Multiplexes together VCRAFT files to give a contiguous, monotonically increasing
    frequency axis
    """

    def __init__(self, vcraft_files, delays=None, default_freq_offset=-1.0):
        self._vcraft_files = sorted(
            vcraft_files, key=lambda fp: (fp.header.card_no, fp.hdr.fpga_id)
        )
        if delays is None:
            delays = {}

        self.file_delays = np.array(
            [
                int(delays.get(os.path.basename(f.fname + ".hdr"), 0))
                for f in vcraft_files
            ]
        )

        freq_str = self.allhdr("FREQS")
        freqs = np.array([map(float, flist.split(",")) for flist in freq_str])

        # there was traditionally a frequency offset of + 1inserted in the header
        # See craft-232
        # undo it if nothing is specified in the header to the contrary

        freq_offset = float(self.hdr_identical("FREQ_OFFSET", default_freq_offset))
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

        self.trigger_frameids = np.array(map(int, self.allhdr("TRIGGER_FRAMEID")))
        self.trigger_mjds = np.array(map(float, self.allhdr("TRIGGER_MJD")))

        self.start_mjds = np.array([f.start_mjd for f in self._files])
        self.start_mjd = max(self.start_mjds)

        self.start_frameids = (
            np.array([f.start_frameid for f in self._files]) - self.file_delays
        )
        self.start_frameid = max(self.start_frameids)
        # BUG! Start MJDs don't account for sample offsets
        self.sample_offsets = self.start_frameid - self.start_frameids

        self.all_overlap_samps = [
            f.nsamps - self.sample_offsets[i] for i, f in enumerate(self._files)
        ]
        self.overlap_nsamps = min(self.all_overlap_samps)

        assert np.all(self.sample_offsets >= 0)
        assert np.all(self.sample_offsets < self.nsamps)
        self.start_mjd = self.start_mjds[np.argmin(self.sample_offsets)]
        assert np.all(
            abs(self.start_mjds - self.trigger_mjds) < 1
        ), "MJD adjustment should be << 1 day"

        beam_ra = np.array(
            map(float, self.allhdr("BEAM_RA", -99))
        ).mean()  # TODO: make craft_vdump write same BEAM_RA for all card/fpgas
        beam_dec = np.array(map(float, self.allhdr("BEAM_DEC", -99))).mean()
        self.beam_pos = (beam_ra, beam_dec)
        pols = np.array([f.pol for f in self._files])
        assert np.all(pols == pols[0])
        self.pol = pols[0]

    @property
    def vcraft_files(self):
        return self._vcraft_files

    @property
    def ant(self):
        return self.check_allhdr("ANT")

    @property
    def antno(self):
        return int(self.check_allhdr("ANTENNA_NO"))

    @property
    def beam(self):
        return int(self.check_allhdr("BEAM"))

    @property
    def nbits(self):
        return int(self.check_allhdr("NBITS"))

    @property
    def mode(self):
        return int(self.check_allhdr("MODE"))

    @property
    def samp_rate(self):
        return float(self.check_allhdr("SAMP_RATE"))

    def get_allhdr(self, cardname, default=None):
        """
        Returns a list of header values for all files, with the given card name
        """
        return [fp.hdr.get(cardname, default) for fp in self.vcraft_files]

    def check_allhdr(self, cardname, default=None):
        """
        Returns the header value for the given cardname.
        Checks the value is the same for all files. If not, it throws a ValueError
        """
        hdr_values = self.get_allhdr(cardname, default)
        hdr_values_unique = set(hdr_values)
        if len(hdr_values_unique) != 1:
            raise ValueError(
                f"Exepected the same header value for {cardname} for all files. Got these values {hdr_values}"
            )

        return hdr_values_unique.pop()

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
        ), "Invalid read request. nsamp={} samp_start ={} differnce{}".format(
            nsamp, samp_start, nsamp - samp_start
        )

        # allocate giant buffer
        d = np.zeros((nsamp, len(self.freqs)), dtype=np.complex64)
        for ifile, f in enumerate(self._files):
            out_chans = self.freqconfig.chanmaps[ifile, :]
            fsamp_start = samp_start + self.sample_offsets[ifile]
            d[:, out_chans] = f.read(fsamp_start, nsamp)
        return d


def mux_by_none(filenames, **kwargs):
    all_files = [VcraftFile(f, **kwargs) for f in filenames]
    all_files.sort(key=lambda f: (f.hdr["ANT"][0], f.hdr["CARD_NO"], f.hdr["FPGA_ID"]))
    return all_files


def mux_by_pol(filenames, delays=None, default_freq_offset=-1.0, **kwargs):
    """
    :return: Dictionary keyened by 'X' or "Y
    """
    all_files = [VcraftFile(f, **kwargs) for f in filenames]
    all_files.sort(key=lambda f: f.pol)
    mux_by_pol = itertools.groupby(all_files, lambda f: f.pol)
    muxes = {
        pol: VcraftMux(list(files), delays, default_freq_offset)
        for pol, files in mux_by_pol
    }

    return muxes


def mux_by_antenna(filenames, delays=None, default_freq_offset=-1.0, **kwargs):
    all_files = [VcraftFile(f, **kwargs) for f in filenames]
    all_files.sort(key=lambda f: f.hdr["ANT"][0])
    ants = [f.hdr["ANT"][0] for f in all_files]
    mux_by_ant = itertools.groupby(all_files, lambda f: f.hdr["ANT"][0])
    muxes = [
        VcraftMux(list(files), delays, default_freq_offset)
        for antname, files in mux_by_ant
    ]
    muxes.sort(key=lambda mux: mux.antno)  # sort by antenna number

    return muxes
