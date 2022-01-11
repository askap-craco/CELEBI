#!/usr/bin/env python
"""
Tied-array beamforming vcraft files, based on "craftcor.py".

Copyright (C) CSIRO 2017
"""

__author__ = [
    "Keith Bannister, CSIRO <keith.bannister@csiro.au>",
    "Hyerin Cho, Curtin University/Swinburne University "
    + "<chyerin1996@gmail.com>",
    "Danica Scott, Curtin University "
    + "<danica.scott@postgrad.curtin.edu.au>",
]

from typing import Tuple, TypeVar

import argparse
import logging
import multiprocessing
import os
import signal
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np
import vcraft
from astropy.coordinates import SkyCoord
from calc11 import ResultsFile
from miriad import MiriadGainSolutions

# constants
C_LIGHT = 299792458.0  # Speed of light (m/s)
CHAN_BWIDTH = 27.0  # Channel bandwidth
OS_NYQ_BWIDTH = 32.0  # Oversampled Nyquist bandwidth
F_OS = OS_NYQ_BWIDTH / CHAN_BWIDTH  # Oversampling factor
NUM_GUARD_CHAN = OS_NYQ_BWIDTH - CHAN_BWIDTH  # Number of guard channels

# generic type
T = TypeVar("T")


# TODO: Cleanup (not in any particular order)
#   IN GENERAL: MAKE THINGS CONSISTENT AND NICE!
#   (1) PEP-008ify everything
#         - Probably do first, easiest to do the others with readable code!
#   (2) Document all functions and classes
#         - Includes docstrings and line comments
#   (3) Ensure only necessary functions remain
#   (4) Remove magic numbers!
#   (5) Make as much as possible compatible with Python 3


class AntennaSource:
    # TODO: (1, 2, 4, 5)
    # Possible refactor: include index and corr as fields

    def __init__(self, vfile):
        """Initialise the AntennaSource object.

        Gets much of the information from the vfile, which is a data
        structure representing the voltages in increasing frequency
        order.

        :param vfile: The vfile for this antenna

        :field vfile: The vfile for this antenna
        :field ant_name: Antenna name. Of the format 'AK**'
        :field antno: Antenna number between 1-36. Note that this is the
                      physical antenna number, and some antennas may not
                      be present.
        :field mjd_start: Start time of the voltages in MJD
        :field trigger_frame: TODO: ?
        :field hdr: Header file for the vfile
        :field all_geom_delays: List of all geometric delays calculated
                                for each channel
        :field all_mjd_mids: List of all central MJDs calculated by the
                             Correlator object
        :field pol: Polarisation of this antenna

        :func do_f_tab: Perform the tied-array beamforming for this
                        antenna.
        """
        # TODO: (5)
        self.vfile = vfile
        self.ant_name = self.vfile.hdr["ANT"][0].lower()
        self.antno = int(self.vfile.hdr["ANTENNA_NO"][0])
        self.mjd_start = self.vfile.start_mjd
        self.trigger_frame = self.vfile.start_frameid
        self.hdr = self.vfile.hdr
        self.all_geom_delays = []
        self.all_mjd_mids = []
        self.pol = self.vfile.pol.lower()
        print(f"antenna {self.ant_name} {self.vfile.freqconfig}")

    def do_f_tab(self, corr: Correlator, i_ant: int) -> np.ndarray:
        """Perform tied-array beamforming and PFB inversion

        :param corr: Correlator object with constants and functions for
            this data set
        :type corr: :class:`Correlator`
        :param i_ant: Index of antenna
        :type i_ant: int
        :return: Beamformed, calibrated, PFB-inverted fine spectrum
        :rtype: np.ndarray
        """
        # TODO: (2, 4, 5)

        self.all_mjd_mids.append(corr.curr_mjd_mid)

        # Initialise output data array
        data_out = np.zeros(
            (corr.n_int, corr.n_fine_chan, corr.n_pol_in), dtype=np.complex64
        )

        n_fine = corr.n_fft - 2 * corr.nguard_chan
        n_samp = corr.n_int * corr.n_fft

        whole_delay, geom_delays_us = self.get_delays(corr, n_samp)

        print(("antenna #: ", i_ant, self.ant_name))
        sample_offset = whole_delay + corr.abs_delay
        frameid = self.vfile.start_frameid + sample_offset
        print(
            "FRAMEID: "
            + str(frameid)
            + ", remainder from 32: "
            + str(frameid % 32)
        )

        """
        To avoid iPFB fractional delay, set FRAMEID such that the 
        remainder is 0
        """

        raw_data = self.vfile.read(sample_offset, n_samp)

        assert raw_data.shape == (
            n_samp,
            corr.ncoarse_chan,
        ), "Unexpected shape from vfile: {} expected ({},{})".format(
            raw_data.shape, n_samp, corr.ncoarse_chan
        )

        for i, chan in enumerate(range(corr.ncoarse_chan)):
            xfguard_f, fine_chan_start, fine_chan_end = process_chan(
                chan, corr, n_fine, raw_data[:, chan], geom_delays_us, i_ant
            )

            """
            Slot xfguard (a trimmed spectrum for this coarse channel) 
            into the corresponding slice of the fine channels
            """
            data_out[:, fine_chan_start:fine_chan_end, 0] = xfguard_f

        return data_out

    def get_delays(
        self, corr: Correlator, n_samp: int
    ) -> "tuple[int, np.ndarray]":
        """Parse and calculate delay for this Antenna from Correlator

        :param corr: Correlator object for this data set
        :type corr: :class:`Correlator`
        :param n_samp: Number of samples in this data set
        :type n_samp: int
        :return: Whole delay and frequency-dependent geometric delay in
            units of microseconds
        :rtype: Tuple[int, :class:`numpy.ndarray`]
        """
        # TODO: (2, 5)

        geom_delay_us, geom_delay_rate_us = corr.get_geom_delay_delayrate_us(
            self
        )

        self.all_geom_delays.append(geom_delay_us)

        fixed_delay_us = corr.get_fixed_delay_usec(self.antno)

        # calculate sample start
        framediff_samp = corr.ref_ant.trigger_frame - self.trigger_frame
        total_delay_samp = framediff_samp
        whole_delay = int(np.round(total_delay_samp))

        # time-dependent geometric delays
        # np.linspace(0, 1, n_samp) == time in units of integrations
        geom_delays_us = (
            geom_delay_us
            + geom_delay_rate_us * np.linspace(0, 1, n_samp)
            - fixed_delay_us
        )

        return whole_delay, geom_delays_us


class FringeRotParams:
    # TODO: (2, 5)

    def __init__(self, corr, ant):
        # TODO: (2, 5)
        mid_data = corr.fringe_rot_data_mid[ant.ant_name]
        self.u = float(mid_data["U (m)"])
        self.v = float(mid_data["V (m)"])
        self.w = float(mid_data["W (m)"])
        self.delay = float(mid_data["DELAY (us)"])
        self.delay_start = float(
            corr.fringe_rot_data_start[ant.ant_name]["DELAY (us)"]
        )
        self.delay_end = float(
            corr.fringe_rot_data_end[ant.ant_name]["DELAY (us)"]
        )
        self.delay_rate = (self.delay_end - self.delay_start) / float(
            corr.n_int
        )
        self.ant = ant
        self.corr = corr

    def __str__(self):
        # TODO: (2, 5)
        s = "FR {} uvw=({},{},{}) m = {} us".format(
            self.ant.ant_name, self.u, self.v, self.w, self.delay
        )
        return s

    __repr__ = __str__


class Correlator:
    """A class to handle processing of raw ASKAP voltages to produce a
    correlated, beamformed fine spectrum via polyphase filterbank (PFB)
    inversion.

    :param ants: A list of AntennaSource objects for the overall dataset
    :type ants: List of AntennaSources
    :param values: Parsed command line arguments
    :type values: :class:`argparse.Namespace`
    :param abs_delay: Absolute delay to be applied uniformly to all
        antennas. Defaults to 0
    :type abs_delay: int
    """

    # TODO: (1, 2, 4, 5)
    def __init__(self, ants, values, abs_delay=0):
        """Constructor function"""
        # TODO: (1, 2, 4, 5)
        self.running = True  # Used to exit gracefully if killed
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)
        self.ants = ants  # All the AntennaSources of the dataset
        self.values = (
            values  # Command line arguments used to initialise this Correlator
        )
        self.pool = None  # If not None, a pool to use for parallelisation
        if self.values.num_threads > 1:
            self.pool = multiprocessing.Pool(processes=values.num_threads)

        self.parse_par_set()

        for ia, a in enumerate(self.ants):
            a.ia = ia  # Set antenna index
            a.antpos = self.get_ant_location(a.antno)  # Get antenna position

        self.abs_delay = abs_delay  # Absolute delay
        self.ref_ant = ants[0]  # Reference Antenna
        self.calc_results = ResultsFile(values.calcfile)  # Parsed calcfile
        self.mjd0 = (
            self.ref_ant.mjd_start
        )  # MJD of data start in reference antenna
        self.frame0 = self.ref_ant.trigger_frame  # ?
        self.n_int = int(values.n_int)  # Number of integrations
        self.n_fft = 64 * values.fft_size  # Number of fine channels to create
        self.nguard_chan = (
            NUM_GUARD_CHAN * values.fft_size
        )  # Number of channels to be chopped off the ends of each coarse channel
        self.oversamp = F_OS  # Oversampling factor
        self.fs = self.oversamp  # samples per microsecond
        self.ncoarse_chan = len(
            self.ref_ant.vfile.freqs
        )  # Number of coarse channels
        self.sideband = -1  # Coarse channel frequency step #?
        self.coarse_chan_bwidth = 1.0  # Coarse channel bandwidth (MHz)
        self.n_fine_per_coarse = (
            self.n_fft - 2 * self.nguard_chan
        )  # Number of fine channels per coarse channel
        self.n_fine_chan = int(
            self.ncoarse_chan * self.n_fine_per_coarse
        )  # Total number of fine channels
        self.fine_chan_bwidth = self.coarse_chan_bwidth / float(
            self.n_fine_per_coarse
        )  # Fine channel width (MHz)
        self.full_bw = (
            self.fine_chan_bwidth * self.n_fine_chan
        )  # Full bandwidth (MHz)
        self.freq_scrunch = values.freq_scrunch  # Frequency averaging factor
        assert self.freq_scrunch >= 1, "Freq scrunch factor must be >= 1"
        assert self.n_fine_per_coarse % self.freq_scrunch == 0, (
            "Freq scrunch must yield an integer number of fine channels "
            "per coarse channel"
        )
        self.n_fine_out_per_coarse = (
            self.n_fine_per_coarse / self.freq_scrunch
        )  # Output number of fine channels per coarse channel
        self.n_fine_out_chan = (
            self.n_fine_out_per_coarse * self.ncoarse_chan
        )  # Output number of fine channels
        self.out_chan_bwidth = self.coarse_chan_bwidth / float(
            self.n_fine_out_per_coarse
        )  # Output fine channel bandwidth
        self.n_pol_in = 1  # Number of polarisations in
        self.f0 = self.ants[0].vfile.freqs[
            0
        ]  # TODO: is this top or bottom of band?
        self.freqs = self.ants[
            0
        ].vfile.freqs  # Coarse channel frequencies (MHz)
        self.freq_mid = (
            self.freqs.mean()
        )  # Central frequency of full band (MHz)
        self.int_time_secs = float(self.n_int * self.n_fft) / (
            self.fs * 1e6
        )  # Integration time in seconds
        self.int_time_days = (
            self.int_time_secs / 86400.0
        )  # Integration time in days
        self.curr_int_no = 0  # Index of currently processing integration
        self.calc_mjd()
        self.get_fringe_rot_data()
        self.pol = self.ants[0].pol  # Polarisation index
        self.parse_mir()

        # TODO: Why are par_set and mir set to None after they're parsed
        self.par_set = None
        self.mir = None

        # Start, mid, and end times of the data in MJD
        self.curr_mjd_start = None
        self.curr_mjd_mid = None
        self.curr_mjd_end = None

        # These are also set to None after being set???
        self.fringe_rot_data_start = None
        self.fringe_rot_data_mid = None
        self.fringe_rot_data_end = None

        logging.debug(
            "F0 %f FINE CHANNEL %f kHz num=%d freqs=%s",
            self.f0,
            self.fine_chan_bwidth * 1e3,
            self.n_fine_chan,
            self.freqs,
        )

    def exit_gracefully(self, signum, frame):
        """Graceful death when e.g. a KeyboardInterrupt happens"""
        self.running = False

    def parse_par_set(self):
        """Parse the delays parset into a dict for later use"""
        self.par_set = {}

        # open the fcm file
        with open(self.values.par_set) as f:
            for line in f:
                if "=" not in line or line.startswith("#"):
                    continue

                name, value = line.strip().split("=")
                name = name.strip()
                value = value.strip()
                self.par_set[name] = value

    def parse_mir(self):
        """Parse calibration solutions and store in mir attribute"""
        self.mir = MiriadGainSolutions(
            self.values.mirsolutions, self.values.aips_c, self.pol, self.freqs
        )

    def get_ant_location(self, antno):
        """Finds the location of the antenna with the given number as
        an ITRF (International Terrestrial Reference Frame) coordinate.

        :param antno: Absolute antenna number (i.e. not index of antenna
            within this dataset)
        :type antno: int
        :return: An ITRF coordinate as a list of three floats [x, y, z]
        :rtype: list
        """
        key = f"common.antenna.ant{antno}.location.itrf"
        value = self.par_set[key]
        location = list(
            map(float, value.replace("[", "").replace("]", "").split(","))
        )
        return location

    def get_fixed_delay_usec(self, antno):
        """Finds the absolute/hardware delay for the antenna with the
        given number in microseconds.

        :param antno: Absolute antenna number (i.e. not index of antenna
            within this dataset)
        :type antno: int
        :return: Absolute delay for this antenna in microseconds
        :rtype: float
        """
        key = f"common.antenna.ant{antno}.delay"
        value = self.par_set[key]
        delay_ns = float(value.replace("ns", ""))
        delay_us = delay_ns / 1e3

        return delay_us

    def get_geom_delay_delayrate_us(self, ant):
        fr1 = FringeRotParams(self, ant)
        fr2 = FringeRotParams(self, self.ref_ant)

        # TODO: There is a discrepancy here, below comment says fr1 is
        #       ref ant, but above suggests fr2 is?
        # fr1: reference antenna
        # Account for effects of Earth's rotation
        delay = fr1.delay_start - fr2.delay_start
        delayrate = fr1.delay_rate - fr2.delay_rate

        return delay, delayrate

    def calc_mjd(self):
        # TODO: (1, 2, 4, 5)
        i = float(self.curr_int_no)
        abs_delay_days = float(self.abs_delay) / 86400.0 / (self.fs * 1e6)
        self.curr_mjd_start = (
            self.mjd0 + self.int_time_days * (i + 0.0) + abs_delay_days
        )
        self.curr_mjd_mid = (
            self.mjd0 + self.int_time_days * (i + 0.5) + abs_delay_days
        )
        self.curr_mjd_end = (
            self.mjd0 + self.int_time_days * (i + 1.0) + abs_delay_days
        )

    def get_calc_results(self, mjd):
        # TODO: (1, 2, 4, 5)
        res = self.calc_results.scans[0].eval_src0_poly(mjd)

        return res

    def get_fringe_rot_data(self):
        # TODO: (1, 2, 4, 5)
        self.fringe_rot_data_start = self.get_calc_results(self.curr_mjd_start)
        self.fringe_rot_data_mid = self.get_calc_results(self.curr_mjd_mid)
        self.fringe_rot_data_end = self.get_calc_results(self.curr_mjd_end)

    def do_tab(self, an=None):
        # TODO: (1, 2, 4, 5)
        # Tied-array beamforming

        n_samp = int(self.n_int)
        n_chan = int(self.ncoarse_chan * self.n_fine_per_coarse)

        sum_aligned = np.zeros(
            (n_samp, n_chan, int(self.n_pol_in)), dtype=np.complex64
        )

        if an == None:  # add all antennas
            print("## Summing up all " + str(len(self.ants)) + " antennas")
            if self.values.tab:
                print("coherent sum")
            else:
                print("incoherent sum")
            for iant, ant in enumerate(self.ants):
                if not self.running:
                    raise KeyboardInterrupt()

                temp = ant.do_f_tab(self, iant)

                if self.values.tab:
                    sum_aligned += temp
                else:
                    sum_aligned += np.power(abs(temp), 2)

            return sum_aligned
        else:
            print("## Operate on only antenna #: " + str(an))
            ant = self.ants[an]
            iant = an
            temp = ant.do_f_tab(self, iant)
            return temp


def _main():
    values = get_args()

    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    antennas = get_antennas(values)

    given_offset = values.offset
    corr = Correlator(antennas, values, abs_delay=given_offset)

    t0 = time.time()
    try:
        print("PERFORMING TIED-ARRAY BEAMFORMING")

        temp = corr.do_tab(values.an)
        fn = values.outfile

        print("saving output to " + fn)

        np.save(fn, temp)
    except Exception as e:
        print("ERROR OCCURRED")
        print(e)
    finally:
        print("craftcor_tab.py running time: " + str(time.time() - t0))
        print("done")


def get_args() -> argparse.Namespace:
    """Parse command line arguments

    :return: Command line argument parameters
    :rtype: argparse.Namespace
    """
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
        default=True,
    )
    parser.add_argument(
        "-o", "--outfile", help="Output fits/.npy file", default="corr.fits"
    )
    parser.add_argument(
        "-c", "--channel", type=int, help="Channel to plot", default=0
    )
    parser.add_argument(
        "-n",
        "--fft-size",
        type=int,
        help="Multiple of 64 channels to make channels - " + "default=1",
        default=1,
    )
    parser.add_argument(
        "-t",
        "--num-threads",
        type=int,
        help="Number of threads to run with",
        default=1,
    )
    parser.add_argument("--calcfile", help="Calc file for fringe rotation")
    parser.add_argument("-w", "--hwfile", help="Hw delay file")
    parser.add_argument("-p", "--par_set", help="Parset for delays")
    parser.add_argument(
        "--show", help="Show plot", action="store_true", default=False
    )
    parser.add_argument(
        "-i",
        "--n_int",
        help="Number of fine spectra to average",
        type=int,
        default=128,
    )
    parser.add_argument(
        "-f",
        "--freq_scrunch",
        help="Frequency average by this factor",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--rfidelay",
        type=int,
        help="Delay in fine samples to add to second component"
        + " to make an RFI data set",
        default=0,
    )
    parser.add_argument(
        "--mirsolutions", help="Root file name for miriad gain solutions"
    )
    parser.add_argument(
        "--aips_c", help="AIPS bandpass polynomial fit coeffs", default=None
    )
    parser.add_argument(
        "--an", type=int, help="Specific antenna", default=None
    )
    parser.add_argument(
        "--offset", type=int, help="FFT offset to add", default=0
    )
    parser.add_argument(dest="files", nargs="+")

    return parser.parse_args()


def process_chan(
    chan: int,
    corr: Correlator,
    n_fine: int,
    chan_raw_data: np.ndarray,
    geom_delays_us: np.ndarray,
    i_ant: int,
) -> "tuple[np.ndarray, int, int]":
    """Process a particular channel of an Antenna.

    Processing includes:
        - Applying fringe rotation to account for the Earth's rotation
        - Trimming the edge of the spectral content of the raw data to
          remove the oversampled regions
        - Applying geometric delays to the data
        - Applying calibration solutions to the data

    :param chan: Channel index
    :type chan: int
    :param corr: Correlator object
    :type corr: :class:`Correlator`
    :param n_fine: Total number of fine channels
    :type n_fine: int
    :param chan_raw_data: Raw, oversampled fine spectrum for this
        channel
    :type chan_raw_data: :class:`numpy.ndarray`
    :param geom_delays_us: Geometric delays as a function of frequency
        in units of microseconds
    :type geom_delays_us: :class:`numpy.ndarray`
    :param i_ant: Antenna index
    :type i_ant: int
    :return: The processed fine spectrum of the channel, and the indexes
        to be used to place the data in the full-band spectrum
    :rtype: Tuple[:class:`np.ndarray`, int, int]
    """

    # TODO: (5)
    # Channel frequency
    centre_freq = corr.freqs[chan]

    """
    Array of fine channel frequencies relative to centre frequency
    Goes from -0.5 to 0.5
    """
    delta_freq = corr.fine_chan_bwidth * (
        np.arange(n_fine, dtype=np.float) - float(n_fine) / 2.0
    )

    # TODO: what is sideband?
    if corr.sideband == -1:
        delta_freq = -delta_freq

    """
    raw_data's shape: (n_samp, corr.ncoarse_chan)
    n_samp = input.i * (64 * input.n)
    """
    x1 = chan_raw_data.reshape(-1, corr.n_fft)

    """
    fixed_delay_us = corr.get_fixed_delay_usec(self.antno)
    fixed_delay_us is contained in fcm.txt for each antenna
    geom_delay_us = corr.get_geometric_delay_delayrate_us(self)[0]
    geom_delay_us accounts for Earth's rotation
    """
    # Fringe rotation for Earth's rotation
    turn_fringe = centre_freq * geom_delays_us
    phasor_fringe = np.exp(2j * np.pi * turn_fringe, dtype=np.complex64)
    x1 *= phasor_fringe

    """
    corr.nguard_chan = NUM_GUARD_CHAN * input.n 
                     = number of fine channels on each side to cut off
    xfguard is xf1 with the ends trimmed off
    """
    xf1 = np.fft.fft(x1, axis=1)
    xf1 = np.fft.fftshift(xf1, axes=1)

    # scale because otherwise it overflows
    xfguard_f = xf1[:, corr.nguard_chan : corr.nguard_chan + n_fine :]

    # Fractional sample phases
    turn_frac = delta_freq * np.mean(geom_delays_us)

    # phasors to rotate the data with constant amplitude = 1
    phasor = np.exp(np.pi * 2j * turn_frac, dtype=np.complex64)

    # get absolute frequencies in gigahertz
    freq_ghz = (centre_freq + delta_freq) / 1e3

    # get the calibration solutions and apply them to the phasors
    mir_cor = corr.mir.get_solution(i_ant, 0, freq_ghz)

    if mir_cor[0] == 0:  # if correction is 0, flag data
        phasor *= 0
    else:
        phasor /= mir_cor

    xfguard_f *= phasor

    # select the channels for this coarse channel
    fine_chan_start = chan * n_fine
    fine_chan_end = (chan + 1) * n_fine

    """
    RECAP
    xfguard is a "dynamic" spectrum of the current coarse 
    channel in only the fine channels not trimmed 
    (oversampled PFB).
    xfguard.shape == (input.i, 64 * input.n)
    """

    return xfguard_f, fine_chan_start, fine_chan_end


def parse_delays(values: argparse.Namespace) -> dict:
    """Parse hardware delays from hwdelays file inferred from imfile

    :param values: Command line argument parameters
    :type values: :class:`argparse.Namespace`
    :return: Dictionary of hardware delays
    :rtype: dict
    """
    # TODO: (2, 4, 5)
    delay_file = values.calcfile.replace(".im", ".hwdelays")

    if os.path.exists(delay_file) is False:
        delay_file = values.hwfile

    delays = {}
    if delay_file is not None and os.path.exists(delay_file):
        with open(delay_file) as dfile:
            for line in dfile:
                bits = line.split()
                if not line.startswith("#") and len(bits) == 2:
                    raw = -int(bits[1])
                    if raw % 8 != 0:  # if it is not a multiple of 8, round
                        new = int(8 * round(float(raw) / 8))
                        print(("hwdelay ", raw, " rounded to ", new))
                        delays[bits[0].strip()] = new
                    else:
                        delays[bits[0].strip()] = raw

        logging.info("Loaded %s delays from %s", len(delays), delay_file)
    else:
        logging.info("No delays loaded. %s does not exist", delay_file)

    return delays


def get_antennas(values: argparse.Namespace) -> "list[AntennaSource]":
    """Create AntennaSource objects for each antenna in data set

    :param values: Command line argument parameters
    :type values: :class:`argparse.Namespace`
    :return: List of all AntennaSource objects for the data set
    :rtype: list[AntennaSource]
    """

    delay_map = parse_delays(values)
    antennas = [
        AntennaSource(mux)
        for mux in vcraft.mux_by_antenna(values.files, delay_map)
    ]

    print(("NUMBER OF ANTENNAS", len(antennas)))

    return antennas


def print_var(name: str, value: T) -> None:
    """Debugging function for output of variable content

    :param name: Name of variable
    :type name: str
    :param value: Variable value
    :type value: T
    """

    base_str = "{} : {}"
    print(base_str.format(name, value))


if __name__ == "__main__":
    _main()
