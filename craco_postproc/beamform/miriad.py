import os
import warnings

import numpy as np
from parse_aips import aipscor
from scipy.interpolate import interp1d


def parse_gpplt(fin):
    # TODO: (1, 2, 3, 4, 5)
    """
    Parse a miriad gpplt exported log file

    :returns: tuple(x, values) where x is an array of strings
    containing x values (dependant on the type of file)
    and values is a (len(x), Nant) numpy array
    """
    with open(fin) as f:

        all_values = []
        curr_values = []
        x = []
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue

            # X value is in the first 16 columns
            # if empty, it's a continuation of previous antennas
            sz = 14
            xfield = line[0:sz].strip()
            bits = list(map(float, line[sz:].split()))
            if xfield == "":
                curr_values.extend(bits)
            else:
                x.append(xfield)
                curr_values = []
                all_values.append(curr_values)
                curr_values.extend(bits)

        v = np.array(all_values)
        # make 2D
        if v.ndim == 1:
            v = v[np.newaxis, :]

        assert v.ndim == 2

    return np.array(x), v


class MiriadGainSolutions:
    # TODO: (1, 2, 3, 4, 5)
    def __init__(self, file_root, bp_c_root=None, pol=None, freqs=None):
        """Loads gpplt exported bandpass and gain calibration solutions.
        Expects 4 files at the following names, produced by miriad gpplt
        with the given options

        Limitations: Currently does no time interpolation - just uses first time
        $file_root.gains.real - yaxis=real
        $file_root.gains.imag - yaxis=imag
        $file_root.bandpass.real - options=bandpass, yaxis=real
        $file_root.bandpass.imag - options=bandpass, yaxis=imag

        """

        # print(file_root)
        # print(bp_c_root)
        if bp_c_root == None and file_root == None:
            print("No bandpass solutions are given")
            # temporary
            g_real = np.full((1, 36), 1, dtype=np.complex64)
            g_imag = np.full((1, 36), 0, dtype=np.complex64)
            self.bp_real = None
        elif bp_c_root == None:
            print("Using MIR bandpass solutions")
            times1, g_real = parse_gpplt(file_root + ".gains.real")
            times2, g_imag = parse_gpplt(file_root + ".gains.imag")
            if 1:  # temp
                g = g_real + 1j * g_imag
                g = 1 / np.conj(g)
                g_real = np.real(g)
                g_imag = np.imag(g)

            assert all(times1 == times2), "Times in gains real/imag dont match"
            assert g_real.shape == g_imag.shape, "Unequal shapes of gain files"

            freqs1, bp_real = parse_gpplt(file_root + ".bandpass.real")
            freqs2, bp_imag = parse_gpplt(file_root + ".bandpass.imag")
            assert all(
                freqs1 == freqs2
            ), "Freqs in bandpass real/imag dont match"

            assert (
                bp_real.shape == bp_imag.shape
            ), "Unequal shapes of bandpass files"
            nant = g_real.shape[1]
            assert (
                bp_real.shape[1] == nant
            ), "Unequal number of antennas in gain and bandpass files"
            self.nant = nant
            self.freqs = np.array(freqs1).astype(np.float)  # convert to float
            self.times = times1  # TODO: Parse times.
            if len(times1) > 1:
                warnings.warn("MiriadSolution can only handle 1 time step")
            self.bp_real = bp_real
            bp_imag *= -1  # -1 for 181112, temp
            self.bp_imag = bp_imag
            self.bp_real_interp = [
                interp1d(
                    self.freqs,
                    bp_real[:, iant],
                    fill_value=(self.bp_real[0, iant], self.bp_real[-1, iant]),
                    bounds_error=False,
                )
                for iant in range(nant)
            ]
            self.bp_imag_interp = [
                interp1d(
                    self.freqs,
                    bp_imag[:, iant],
                    fill_value=(self.bp_imag[0, iant], self.bp_imag[-1, iant]),
                    bounds_error=False,
                )
                for iant in range(nant)
            ]
            self.bp_coeff = None
        else:
            print("Using AIPS bandpass solutions")
            if "polyfit_coeff" in bp_c_root:  # AIPS polyfit coefficient
                self.bp_coeff = np.load(bp_c_root)
            else:

                nant = None
                nfreq = None
                with open(bp_c_root) as fl:
                    for line in fl:
                        if "NAXIS2" in line:
                            nant = int(line.split()[2])
                            if "190608" in bp_c_root:
                                nant -= 2  # exclude ak31 and ak32 from calibration solutions (from Adam)
                        if "TFDIM11" in line:
                            nfreq = int(line.split()[2])
                if nant is None or nfreq is None:
                    print(
                        "WARNING! nant or nfreq not assigned while parsing AIPS bandpass"
                    )
                fmax = freqs[0] + 0.5  # in MHz
                bw = len(freqs)  # in MHz
                self.freqs = (
                    -np.arange(float(nfreq)) / nfreq * bw
                    + fmax
                    - float(bw) / nfreq / 2
                ) / 1e3  # reassign freqs in GHz
                self.bp_real = np.full(
                    (nfreq, nant), np.nan, dtype=np.complex64
                )
                self.bp_imag = np.full(
                    (nfreq, nant), np.nan, dtype=np.complex64
                )
                g_real = np.full((1, nant), np.nan, dtype=np.complex64)
                g_imag = np.full((1, nant), np.nan, dtype=np.complex64)

                # look for a README and get fring, selfcal filenames
                drcal = os.path.dirname(bp_c_root)
                import glob

                readme = glob.glob(drcal + "/README*")
                if len(readme) == 1:
                    with open(readme[0]) as fl:
                        for line in fl:
                            if "delays" in line and ".sn.txt" in line:
                                fring_f = drcal + "/" + line.split()[0]
                            if "selfcal" in line and ".sn.txt" in line:
                                sc_f = drcal + "/" + line.split()[0]
                                print(f"FOUND SELFCAL FILE: {sc_f}")
                else:
                    print("No or multiple readme file exists for AIPS")
                    fring_f = bp_c_root.replace(
                        bp_c_root.split("/")[-1], "delays.sn.txt"
                    )  # ("bandpasses.bp.txt","delays.sn.txt")
                    sc_f = bp_c_root.replace(
                        bp_c_root.split("/")[-1], "selfcal.sn.txt"
                    )
                    # ("bandpasses.bp.txt","selfcal.sn.txt")

                specific_frb = None
                if "190608" in bp_c_root:
                    specific_frb = "190608"
                aips_cor = aipscor(
                    fring_f, sc_f, bp_c_root, specific_frb=specific_frb
                )
                for iant in range(nant):
                    bp = aips_cor.get_phase_bandpass(iant, pol)
                    bp = np.fliplr([bp])[0]  # decreasing order

                    # fring delay
                    delta_t_fring_ns = (
                        aips_cor.get_delay_fring(iant, pol) * 1e9
                    )
                    phases = delta_t_fring_ns * self.freqs
                    phases -= phases[
                        int(len(phases) / 2)
                    ]  # TODO! READ THE REFERENCE FREQUENCY AND SET TO THAT REFERENCE
                    # print(np.shape(phases),np.shape(bp))
                    bp *= np.exp(np.pi * 2j * phases, dtype=np.complex64)
                    if 0:  # temp
                        print(
                            "Polarization fring file exists, also applying this delay"
                        )
                        if pol == "y":
                            phase_offset = (
                                -9.79084693013
                            )  # 6.66130639185 #10.2025399576 #10.347862410779019#(offsetlarger)# temp, got from VELA
                            polfring_delay = 1.10488043784e-9  # -1.10632473155e-9 #-1.13974507853e-9#-1.1568853748356676e-9#(delaylarger)
                        else:
                            phase_offset = 0
                            polfring_delay = 0
                        print((polfring_delay * 1e9, phase_offset))  # temp
                        polphases = (
                            2 * np.pi * polfring_delay * 1e9 * self.freqs
                            + phase_offset
                        )
                        bp *= np.exp(polphases * 1j, dtype=np.complex64)
                    try:
                        g = aips_cor.get_phase_fring(
                            iant, pol
                        ) * aips_cor.get_phase_selfcal(iant, pol)
                        g = 1 / g  # inverse of gain
                    except Exception as e:
                        print(e)
                        g = 0
                    # g = np.conj(g)
                    bp = np.conj(bp)
                    if "180924" in bp_c_root and iant == 15:
                        print(
                            "FRB 180924 antenna 15 for AIPS signal negated. fixing this..."
                        )
                        bp *= -1
                    if "190102" in bp_c_root:
                        print(
                            "FRB 190102 has bandpass magnitude inversed. fixing this..."
                        )
                        bp = 1 / np.conj(bp)
                    self.bp_real[:, iant] = np.real(bp)
                    self.bp_imag[:, iant] = np.imag(bp)
                    g_real[0, iant] = np.real(g)
                    g_imag[0, iant] = np.imag(g)
                # self.freqs = (-np.arange(2016)/2016.*336+1488.5-1/12.)/1e3 # reassign freqs in GHz
                # self.bp_real, self.bp_imag = parse_aips_bp(bp_c_root, pol)
                self.bp_real_interp = [
                    interp1d(
                        self.freqs,
                        self.bp_real[:, iant],
                        fill_value=(
                            self.bp_real[0, iant],
                            self.bp_real[-1, iant],
                        ),
                        bounds_error=False,
                    )
                    for iant in range(nant)
                ]
                self.bp_imag_interp = [
                    interp1d(
                        self.freqs,
                        self.bp_imag[:, iant],
                        fill_value=(
                            self.bp_imag[0, iant],
                            self.bp_imag[-1, iant],
                        ),
                        bounds_error=False,
                    )
                    for iant in range(nant)
                ]
                self.bp_coeff = None

        self.g_real = g_real
        self.g_imag = g_imag
        # arrays indexed by antenna index

    def get_solution(self, iant, time, freq_ghz):
        # TODO: (1, 2, 3, 4, 5)
        """
        Get solution including time and bandpass
        i_ant - antenna index
        time - some version of time. Ignored for now
        freq_ghz - frequency float in Ghz
        """
        if (
            self.bp_real is None
        ):  # this means no bandpass/gain solution was passed
            bp_value = np.array([1])
        elif self.bp_coeff is not None:  # Use AIPS polyfit coefficient
            bp_fit = np.poly1d(self.bp_coeff[iant, 0, :]) + 1j * np.poly1d(
                self.bp_coeff[iant, 1, :]
            )
            bp_value = bp_fit(freq_ghz * 1e3)
        else:  # AIPS polyfit coefficient doesn't exist. Use Miriad/AIPS bandpass interpolation
            f_real = self.bp_real_interp[iant](freq_ghz)
            f_imag = self.bp_imag_interp[iant](freq_ghz)
            bp_value = f_real + 1j * f_imag

        g_value = self.g_real[0, iant] + 1j * self.g_imag[0, iant]
        total_value = bp_value * g_value

        return total_value

    def plot(self):
        # TODO: (1, 2, 3, 4, 5)
        """fig, ax = pylab.subplots(3,3)
        ax = ax.flatten()
        for i in xrange(min(9, self.nant)):
            freq_ghz = self.freqs
            sol = self.get_solution(i, 0, freq_ghz)
            ax[i].plot(freq_ghz, abs(sol))
            ax2 = ax[i].twinx()
            ax2.plot(freq_ghz, np.degrees(np.angle(sol)), 'ro')
            ax[i].set_xlabel('Freq(GHz)')
        """
        return  # added by hyerin
