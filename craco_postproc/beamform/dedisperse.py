#!/usr/bin/env python3

import time

import numpy as np


def _main():
    start = time.time()
    args = get_args()
    f = load(args.f)
    f_dd = dedisperse(f, args.DM, args.f0, args.bw)
    save(f_dd, args.o)
    end = time.time()
    print(f"dedisperse.py finished in {end-start} s")


def get_args():
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

    parser = ArgumentParser(
        description="Coherently dedisperses given fine spectrum",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", help="Fine spectrum file to dedisperse")
    parser.add_argument(
        "--DM",
        type=float,
        help="Dispersion measure to dedisperse to in pc/cm3",
    )
    parser.add_argument(
        "--f0", type=float, help="Central frequency of observation in MHz"
    )
    parser.add_argument(
        "--bw", type=float, help="Bandwidth of observation in MHz", default=336
    )
    parser.add_argument(
        "-o", help="Output file to save dedispersed spectrum to"
    )
    return parser.parse_args()


def load(fname):
    print(f"Loading {fname}")
    return np.load(fname)


def dedisperse(f, DM, f0, bw):
    """
    Takes heavy inspiration from Hyerin Cho's coh_dedisp function from freq2time.py
    """
    # print('Dedispersing')

    n_sam = len(f)

    """
	This value of k_DM is not the most precise available. It is used because to alter the commonly-used
	value would make pulsar timing very difficult. Also, to quote Hobbs, Edwards, and Manchester 2006:
		...ions and magnetic fields introduce a rather uncertain correction of the order of a part in
		10^5 (Spitzer 1962), comparable to the uncertainty in some measured DM values...
	"""
    k_DM = 2.41e-4

    f_min = f0 - float(bw) / 2
    f_max = f0 + float(bw) / 2

    # TODO: figure out why Hyerin put f_max and f_min in this order
    freqs = np.linspace(f_max, f_min, n_sam)

    dedisp_phases = np.exp(
        2j * np.pi * DM / k_DM * ((freqs - f0) ** 2 / f0 ** 2 / freqs * 1e6)
    )

    f *= dedisp_phases

    return f


def save(f_dd, fname):
    print(f"Saving {fname}")
    np.save(fname, f_dd)


if __name__ == "__main__":
    _main()
