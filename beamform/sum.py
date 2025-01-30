"""
Sums time series together for a given FRB

THIS ENTIRE FILE IS PROBABLY REDUNDANT - CAN PROCESS ALL ANTENNAS AND
SUM THEM IN craftcor_tab.py
"""

import glob
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np


def _main():
    start = time.time()
    args = get_args()

    fnames = find_files(args.f_dir, args.f, args.p)

    sum = do_sum(fnames)
    save(sum, args.o)
    end = time.time()
    print(f"sum.py completed in {end-start} s")


def get_args():
    parser = ArgumentParser(
        description="Sums fine channel spectra together for a given FRB",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--f_dir",
        type=str,
        help="Directory where the spectra files are stored",
    )
    parser.add_argument("-f", type=str, help="FRB name")
    parser.add_argument("-p", type=str, help="Polarisation to sum")
    parser.add_argument("-o", type=str, help="Output filename")

    return parser.parse_args()


def find_files(f_dir, FRB, p):
    glob_str = f"{f_dir}/{FRB}_*_{p}_f_filtered.npy"
    return glob.glob(glob_str)



def do_sum(fnames):
    # Sum one file at a time, since they're very large
    # Numpy's mmap_mode helps by allowing us to do this without moving them into memory

    print(fnames)

    # Initialise summation array with the shape of the first file
    print("Initialising sum array...")
    f = np.load(fnames[0])
    print(f"f shape: {f.shape}")
    sum_arr = np.zeros(f.shape, dtype=f.dtype)

    n_antennas = 0

    # Go over the rest of the files, summing their contents into sum_arr
    for fname in fnames:
        print(fname)
        # mmap_mode='r' keeps the loaded array on disk, instead of loading it into memory
        new_arr = np.load(fname, mmap_mode="r")
        print(f"Array Size: {new_arr.shape}")
        if not np.isnan(new_arr).any():
            sum_arr += new_arr
            n_antennas += 1

    print(f"Number of good antennas: {n_antennas}")

    return sum_arr


def save(arr, outfile):
    print("Saving...")
    np.save(outfile, arr)


if __name__ == "__main__":
    _main()
