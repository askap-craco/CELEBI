#!/usr/bin/env python
"""
reconstructs FFFF via de-rippling, coherent de-dispersing, and IFFT-ing.
"""

import os

import numpy as np
import scipy.fftpack as fft
from scipy.interpolate import interp1d


def deripple(FFFF, fftLength=1048576, quiet=False, bw=336):
    if not quiet:
        print("derippling....")
    FFFF = FFFF[0, :, 0]

    # ASKAP Parameters
    N = 1536
    OS_De = 27.0
    OS_Nu = 32.0
    passbandLength = int(((fftLength / 2) * OS_De) / OS_Nu)

    # de-ripple coefficients
    dr_c_file = "../Calibration/deripple_res6_nfft" + str(fftLength) + ".npy"
    if os.path.exists(dr_c_file) == False:
        from generate_deripple import generate_deripple

        print("No derippling coefficient found. Generating one...")
        generate_deripple(fftLength, 6)
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


def coh_dedisp(FFFF, DM, f_mid=1320.5, bw=336, quiet=False):
    nSam = len(FFFF)
    k_DM = 4148.808

    # ASKAP Parameters
    f_start = f_mid - float(bw) / 2  # 1153.
    f_stop = f_mid + float(bw) / 2  # 1488.

    if not quiet:
        print("dedispersing....")
    dedisperse = np.exp(
        2j
        * np.pi
        * DM
        * k_DM
        * np.array(
            [
                (f - f_mid) ** 2 / f_mid ** 2 / f * 1e6
                for f in np.linspace(f_stop, f_start, nSam)
            ]
        )
    )
    # print('dedispersing wrong....')
    # dedisperse = np.exp(2j*np.pi*DM*4150*np.array([(1/f**2-1/f_mid**2)*f*1e6 for f in np.linspace(f_stop,f_start,nSam)]))

    FFFF *= dedisperse

    return FFFF


def ifft_long(FFFF, quiet=False):
    if not quiet:
        print("ifft-ing....")
    t_series = fft.ifft(fft.fftshift(FFFF))
    return t_series


def reconstruct(fn, fftLength, DM, f0=1320.5, bw=336, quiet=False):
    FFFF = np.load(fn)
    FFFF = deripple(FFFF, fftLength, quiet, bw)
    FFFF = coh_dedisp(FFFF, DM, f0, bw, quiet)
    return FFFF


def _main():
    import time
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

    t0 = time.time()
    parser = ArgumentParser(
        description="Script description",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--fn", help="FFFF file directory", default=None)
    parser.add_argument(
        "-d", "--DM", type=float, help="Dispersion measure", default=None
    )
    parser.add_argument(
        "--f0", type=float, help="Central frequency", default=1320.5
    )
    parser.add_argument(
        "-o", "--outfile", help="Output time series file directory"
    )
    parser.add_argument("--bw", type=int, help="Bandwidth", default=336)
    parser.add_argument(
        "-l", "--fftlength", type=int, help="FFT length", default=1048576
    )
    parser.add_argument(
        "--no_dr", help="Don't deripple", action="store_true", default=False
    )
    parser.add_argument(
        "--no_dd", help="Don't dedisperse", action="store_true", default=False
    )
    parser.add_argument(
        "--no_ifft", help="Don't ifft", action="store_true", default=False
    )
    parser.add_argument(
        "-q", help="Quiet mode", action="store_true", default=False
    )
    values = parser.parse_args()

    if values.no_dr or values.no_dd or values.no_ifft:
        FFFF = np.load(values.fn)
        if values.no_dr and values.no_dd:
            print("No derippling, no dedispersing.")
            FFFF = FFFF[0, :, 0]
        elif values.no_dr:
            print("No derippling")
            FFFF = FFFF[0, :, 0]
            FFFF = coh_dedisp(
                FFFF, values.DM, f_mid=values.f0, bw=values.bw, quiet=values.q
            )
        else:
            print("No dedispersing")
            FFFF = deripple(
                FFFF, values.fftlength, quiet=values.q, bw=values.bw
            )
        if values.no_ifft:
            print("No ifft. Saving FFFF")
            t_series = FFFF
        else:
            t_series = ifft_long(FFFF, quiet=values.q)
    else:
        print("Reconstructing...")
        FFFF_dd = reconstruct(
            values.fn,
            values.fftlength,
            values.DM,
            f0=values.f0,
            quiet=values.q,
            bw=values.bw,
        )

    if values.outfile is not None:
        print("output saved to: " + values.outfile)
        np.save(values.outfile, FFFF_dd)
    print("freq2time.py running time: " + str(time.time() - t0))


if __name__ == "__main__":
    _main()
