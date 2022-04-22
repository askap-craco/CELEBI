import sys
import time

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
            calculate_stokes(args, x_ds, y_ds, args.o, "dynspec")

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


def calculate_stokes(args, x, y, outfile, type):
    # lambda functions for each of the Stokes parameters
    stokes = {
        "i": lambda x, y: np.abs(x) ** 2 + np.abs(y) ** 2,
        "q": lambda x, y: np.abs(x) ** 2 - np.abs(y) ** 2,
        "u": lambda x, y: 2 * np.real(np.conj(x) * y),
        "v": lambda x, y: 2 * np.imag(np.conj(x) * y),
    }

    stk_args = [args.I, args.Q, args.U, args.V]

    for idx, stk in enumerate(["i", "q", "u", "v"]):
        if stk_args[idx]:
            print(f"Calculating {stk} {type}")
            par = stokes[stk](x, y)
            if type == "dynspec":
                par_norm = normalise(par)
                del par
                par = par_norm.transpose()
                del par_norm

            save(par, outfile, stk, type)


def normalise(ds):
    """
    Normalises a dynamic spectrum along frequency channels to cancel out RFI

    :param ds: Input dynspec to normalise
    """
    T, F = ds.shape
    out_ds = ds.copy()
    norm = lambda x: (x - np.mean(x)) / np.std(x)
    for f in range(F):
        out_ds[:, f] = norm(ds[:, f])
    del ds
    return out_ds


if __name__ == "__main__":
    _main()
