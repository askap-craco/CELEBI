#!/usr/bin/env python
"""
Header utility classes

Copyright (C) CSIRO 2016
"""
__author__ = "Keith Bannister <keith.bannister@csiro.au>"

import warnings

import numpy as np


class FreqConfig:
    """Containiner class containing information about the frequency config for each
    beamformer chassis of an antenna
    """

    def __init__(self, freqs, reverse=True):
        """Makes a freqconfig  configuration that unifies the beamformer frequencies (which are non-continuous)
        Into a start frequency, bandwidth and a set of channel aand frequency maps.
        Some effort is made to handle gaps

        :freqs: an (Nbeamforms, NchansPerBeamformer) array of frequencies in MHz
        :reverse: if True, make the 0th channel the highest one, and the bandwidth negative, Otherwise the 0th channel
        is the lowest one, and the bandwidth is positive
        """

        nbf, nfreqsperbf = freqs.shape
        assert nbf >= 1

        allfreqs = np.sort(freqs.flat)  # also flattens array
        if reverse:
            allfreqs = allfreqs[::-1]

        assert len(np.unique(allfreqs)) == len(
            allfreqs
        ), f"Duplicate frequencies = {str(allfreqs)}"
        freq_diffs = allfreqs[1:] - allfreqs[0:-1]
        bw = np.median(freq_diffs)
        startfreq = allfreqs[0]
        nchan = len(allfreqs)

        if reverse:
            assert bw < 0, "Expected negative bw. Got: %s" % bw
            assert startfreq == max(freqs.flat)
        else:
            assert bw > 0, "Expected positive bw. Got: %s" % bw
            assert startfreq == min(freqs.flat)

        assert bw != 0

        self.nchan = nchan
        self.bw = bw
        self.freq = startfreq
        self.chanmaps = np.empty((nbf, nfreqsperbf), dtype=np.int32)
        self.freqmaps = np.empty((nbf, nfreqsperbf), dtype=np.float32)

        for ibf in range(nbf):
            bf_freq = freqs[ibf, :]
            bfsort = np.sort(bf_freq)
            if reverse:
                bfsort = bfsort[::-1]

            chans = np.round((bf_freq - startfreq) / bw)
            assert np.all(chans >= 0)

            self.chanmaps[ibf, :] = chans
            self.freqmaps[ibf, :] = bfsort

        # check we haven't got duplicate channels and/or frequencies
        assert (
            len(np.unique(self.chanmaps)) == nchan
        ), "Duplicate channel in chanmap"
        assert (
            len(np.unique(self.freqmaps)) == nchan
        ), "Duplicte freq in freqmap"

        self.span_chan = self.chanmaps.max() - self.chanmaps.min()
        self.span_freq = max(allfreqs) - min(allfreqs)

        self.nchan_span = self.span_chan + 1

        assert (
            self.nchan_span >= nchan
        ), f"Chan span {self.span_chan} must be > nchan {nchan}"

        self.has_gaps = self.span_chan > self.nchan
        self.freqs = freqs

    def __str__(self):
        s = "FREQCONFIG freq {} bw {} nchan {} span_chan {} span_freq {} nchan_span {} gaps? {}".format(
            self.freq,
            self.bw,
            self.nchan,
            self.span_chan,
            self.span_freq,
            self.nchan_span,
            self.has_gaps,
        )
        return s

    __repr__ = __str__

    @staticmethod
    def load_from_dada_header(dada_header, cards=None, freqname="FREQS"):
        from crafthdr import DadaHeader

        h = dada_header
        nbf = int(h.get_value("NUM_BEAMFORMERS"))
        freqs = []
        card_data = {}
        for i in range(nbf):
            addr = h.get_value(f"BEAMFORMER{i}_ADDR")
            cardno = int(h.get_value(f"BEAMFORMER{i}_CARDNO"))
            myfreqs = list(
                map(float, h.get_value(f"BEAMFORMER{i}_{freqname}").split(","))
            )

            if len(myfreqs) == 0:
                raise ValueError(f"no frequencies for beamformerer ibf {i}")

            card_data[cardno] = (addr, cardno, myfreqs)
            freqs.append(myfreqs)

        if cards is not None:
            freqs = []
            # Make order important
            for card in cards:
                freqs.append(card_data[card][2])

        freqs = np.array(freqs)
        freqconfig = FreqConfig(freqs)
        return freqconfig

    @staticmethod
    def load_from_pv(ant, cards):
        from epics import PV

        all_freqs = []

        for card in cards:
            chanpv = PV(f"{ant}:abf:c{card:02d}:F_setChanFreq:chans")
            pvchans = chanpv.get()

            # Poke it to get it to read
            pvname = f"{ant}:abf:c{card:02d}:fpgaSubbands"
            subbandpv_proc = PV(pvname + ".PROC")
            ret = subbandpv_proc.put(1)
            assert ret == 1, "Invalid return value for PV %s" % subbandpv_proc

            # Now read that puppy
            subbandpv = PV(pvname)
            for i in range(10):
                subbands = subbandpv.get()
                if len(subbands) > 0:
                    break
                print("Retry", pvname)

            # need to do this because the other bits arent part of the actual subband.
            subbands = subbands & 0x7F

            # Numberof fpgas
            assert (
                len(subbands) == 6
            ), f"Invalid subbands: {subbands} for {subbandpv}"
            assert np.all(
                subbands >= 0
            ), f"Invalid subbands: {subbands} for {subbandpv}"
            pvfreqs = []
            for sbid in subbands:
                for f in range(8):
                    channo = sbid * 8 + f
                    pvfreqs.append(pvchans[channo])

            all_freqs.append(pvfreqs)

        all_freqs = np.array(all_freqs)

        warnings.warn(
            "Adding a henious hack found by looking at pulsar data where we appear to be off by 1 MHz"
        )
        all_freqs += 1
        freqconfig = FreqConfig(all_freqs)

        return freqconfig
