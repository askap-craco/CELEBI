#!/usr/bin/env python

import os
import time

import numpy as np
from scipy.interpolate import interp1d


def _main():
    start = time.time()
    args = get_args()
    sum_f = load(args.f)
    sum_f_deripple = deripple(sum_f, args.coeffs, fftLength=args.l)
    save(sum_f_deripple, args.o)
    end = time.time()
    print(f"deripple.py finished in {end-start} s")


def get_args():
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

    parser = ArgumentParser(
        description="Deripples fine spectrum",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", help="Fine spectrum file")
    parser.add_argument("-l", type=int, help="FFT length")
    parser.add_argument("-o", help="Output file")
    parser.add_argument(
        "-c", "--coeffs", help="Directory for derippling coefficients"
    )
    return parser.parse_args()


def load(fname):
    return np.load(fname)


def deripple(FFFF, coeff_dir, fftLength=1048576, bw=336):
    """
    deripple input fine spectrum
    Taken from freq2time.py, written by Hyerin Cho

    :param FFFF: Input fine spectrum
    :param fftLength: TODO: Figure out what this exactly is for
    :param bw: Bandwidth of observation in MHz
    :return: Derippled fine spectrum
    """
    print("derippling....")
    FFFF = FFFF[0, :, 0]

    # ASKAP Parameters
    N = 1536
    OS_De = 27.0
    OS_Nu = 32.0
    passbandLength = int(((fftLength / 2) * OS_De) / OS_Nu)

    # de-ripple coefficients
    dr_c_file = coeff_dir + "/deripple_res6_nfft" + str(fftLength) + ".npy"
    if os.path.exists(dr_c_file) == False:
        from generate_deripple import generate_deripple

        print("No derippling coefficient found. Generating one...")
        generate_deripple(fftLength, 6, coeff_dir)
    print(f"loading {dr_c_file}")
    temp = np.load(dr_c_file)
    print(f"{dr_c_file} loaded!")
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


def save(f, fname):
    np.save(fname, f)


if __name__ == "__main__":
    _main()
