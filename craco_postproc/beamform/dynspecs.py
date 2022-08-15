import sys
import time

import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np


def _main():
    start = time.time()

    args = get_args()

    x = load(args.x)
    y = load(args.y)

    if not (args.ds or args.t):
        print("Not generating dynamic spectra or time series - exiting!")
        sys.exit(1)

    if not (args.X or args.Y or args.I or args.Q or args.U or args.V):
        print("Not saving any polarisation or Stokes data - exiting!")
        sys.exit(2)

    if args.t:
        if args.I or args.Q or args.U or args.V:
            print("Calculating Stoke parameters")
            calculate_stokes(args, x, y, args.o, "t")

    if args.ds:
        print("Generating x and y dynamic spectra")
        x_ds = generate_dynspec(x)
        y_ds = generate_dynspec(y)

        if args.X:
            save(x_ds, args.o, "x", "dynspec")

        if args.Y:
            save(y_ds, args.o, "y", "dynspec")

        if args.I or args.Q or args.U or args.V:
            print("Calculating Stoke parameters")
            stokes_fnames = calculate_stokes(args, x_ds, y_ds, args.o, "dynspec")
            with open("dynspec_fnames.txt", "w") as f:
                for fname in stokes_fnames:
                    f.write(f"{fname}\n")
        
    end = time.time()
    print(f"dynspecs.py finished in {end - start} s")


def get_args():
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

    parser = ArgumentParser(
        description="Calculates Stokes parameters and generates dynamic spectra"
    )
    parser.add_argument("-x", help="x polarisation time series")
    parser.add_argument("-y", help="y polarisation time series")
    parser.add_argument(
        "-o",
        help="Output file name. Should have a ! that will be replaced with {x, y, i, q, u, v} and a @ that will be replaced with {t, dynspec}",
    )
    parser.add_argument(
        "-ds",
        action="store_true",
        default=False,
        help="Generate dynamic spectra",
    )
    parser.add_argument(
        "-t", action="store_true", default=False, help="Generate time series"
    )
    parser.add_argument(
        "-X",
        action="store_true",
        default=False,
        help="Save X polarisation data",
    )
    parser.add_argument(
        "-Y",
        action="store_true",
        default=False,
        help="Save Y polarisation data",
    )
    parser.add_argument(
        "-I", action="store_true", default=False, help="Save Stokes I data"
    )
    parser.add_argument(
        "-Q", action="store_true", default=False, help="Save Stokes Q data"
    )
    parser.add_argument(
        "-U", action="store_true", default=False, help="Save Stokes U data"
    )
    parser.add_argument(
        "-V", action="store_true", default=False, help="Save Stokes V data"
    )
    return parser.parse_args()


def load(fname):
    print(f"Loading {fname}")
    return np.load(fname)


def generate_dynspec(t_ser):
    """
    Creates a dynamic spectrum at the highest time resolution from the given time series.
    Assumes 336 frequency channels.

    :param t_ser: input time series of voltages
    :return: dynamic spectrum of voltages
    """
    n = 336
    dynspec = np.zeros((int(t_ser.shape[0] / n), n), dtype=np.complex64)
    for i in range(int(t_ser.shape[0] / n)):
        dynspec[i, :] = np.fft.fft(t_ser[i * n : (i + 1) * n])

    return dynspec


def save(arr, fname, id, type):
    save_fname = fname.replace("!", id).replace("@", type)
    print(f"Saving {save_fname}")
    np.save(save_fname, arr)
    return save_fname


def calculate_stokes(args, x, y, outfile, type):
    # lambda functions for each of the Stokes parameters
    stokes = {
        "I": lambda x, y: np.abs(x) ** 2 + np.abs(y) ** 2,
        "Q": lambda x, y: np.abs(x) ** 2 - np.abs(y) ** 2,
        "U": lambda x, y: 2 * np.real(np.conj(x) * y),
        "V": lambda x, y: 2 * np.imag(np.conj(x) * y),
    }

    stk_args = [args.I, args.Q, args.U, args.V]
    fnames = []

    pars = ["I"] if type == "t" else ["I", "Q", "U", "V"]

    for idx, stk in enumerate(pars):
        if stk_args[idx]:
            print(f"Calculating {stk} {type}")
            par = stokes[stk](x, y)
            if type == "dynspec":
                this_means, this_stds = get_norm(par)
                if stk == "I":
                    stds = this_stds

                par_norm = (par - this_means)/stds
                del par
                par = par_norm.transpose()
                del par_norm

            fnames.append(save(par, outfile, stk, type))
    
    return fnames


def get_norm(ds):
    """
    Gets normalisation parameters to apply to dynamic spectra to
    normalise them

    :param ds: Input dynspec to normalise
    """
    T, F = ds.shape
    means = np.tile(np.mean(ds, axis=0), [T, 1])
    stds = np.tile(np.std(ds, axis=0), [T, 1])
    return means, stds


if __name__ == "__main__":
    _main()
