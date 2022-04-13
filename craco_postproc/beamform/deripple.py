#!/usr/bin/env python

import os
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

import numpy as np
from scipy.interpolate import interp1d
from scipy import fft
from scipy import io


def _main():
    start = time.time()
    args = get_args()
    sum_f = np.load(args.f)
    sum_f_deripple = deripple(sum_f, args.coeffs, args.l, args.bw)
    np.save(args.o, sum_f_deripple)
    end = time.time()
    print(f"deripple.py finished in {end-start} s")


def get_args() -> Namespace:
    """Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description="Deripples fine spectrum",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", help="Fine spectrum file")
    parser.add_argument("-l", type=int, help="FFT length")
    parser.add_argument("-o", help="Output file")
    parser.add_argument(
        "-c", "--coeffs", help="Directory to store derippling coefficients"
    )
    parser.add_argument(
        "--bw", type=float, help="Spectrum bandwidth in MHz", default=336
    )
    return parser.parse_args()


def deripple(
    FFFF: np.ndarray, coeff_dir: str, fftLength: int, bw: float
) -> np.ndarray:
    """Deripple a fine spectrum.

    :param FFFF: Fine spectrum to be derippled
    :type FFFF: :class:`np.ndarray`
    :param coeff_dir: Directory where derippling coefficients are stored
    :type coeff_dir: str
    :param fftLength: Probably can be inferred from shape of FFFF
    :type fftLength: int
    :param bw: Bandwidth of spectrum in MHz
    :type bw: float
    :return: Derippled fine spectrum
    :rtype: :class:`np.ndarray`
    """
    print("derippling....")
    FFFF = FFFF[0, :, 0]

    # ASKAP Parameters
    N = 1536
    OS_De = 27.0
    OS_Nu = 32.0
    passbandLength = int(((fftLength / 2) * OS_De) / OS_Nu)

    # de-ripple coefficients
    dr_c_file = f"{coeff_dir}/deripple_res6_nfft{fftLength}.npy"
    if os.path.exists(dr_c_file) == False:
        print("No derippling coefficient found. Generating one...")
        generate_deripple(fftLength, 6, coeff_dir)

    temp = np.load(dr_c_file)
    print("Interpolating...")

    interp = interp1d(6 * np.arange(len(temp)), temp)
    print("Calculating deripple...")

    deripple = np.ones(passbandLength + 1) / abs(
        interp(np.arange(passbandLength + 1))
    )

    for chan in range(bw):
        print(chan)
        for ii in range(passbandLength):
            FFFF[ii + chan * passbandLength * 2] = (
                FFFF[ii + chan * passbandLength * 2]
                * deripple[passbandLength - ii]
            )
            FFFF[passbandLength + ii + chan * passbandLength * 2] = (
                FFFF[passbandLength + ii + chan * passbandLength * 2]
                * deripple[ii]
            )

    return FFFF


def generate_deripple(nfft, res, dir):
    N = 1536
    OS_De = 27.0
    OS_Nu = 32.0
    dr = dir + "/ADE_R6_OSFIR.mat"
    h = io.loadmat(dr)["c"][0]
    passbandLength = int(((nfft / 2) * OS_De) / OS_Nu)
    multiple = int(1536 / res)

    # Modifications: pre-pad h with zeros
    #                use scipy's fft instead of numpy's, as np's was segfaulting for huuuuge multiple*passbandLength*2
    h_0 = np.zeros(multiple * passbandLength * 2)
    h_0[: h.shape[0]] = h
    temp = abs(fft.fft(h_0))  # ,multiple*passbandLength*2))
    print((
        "saving {}".format(
            dir + "/deripple_res" + str(res) + "_nfft" + str(nfft)
        )
    ))
    np.save(dir + "/deripple_res" + str(res) + "_nfft" + str(nfft), temp)


if __name__ == "__main__":
    _main()
