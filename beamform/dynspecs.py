# from mmap import mmap
import sys
import time

# import matplotlib.gridspec as gs
# import matplotlib.pyplot as plt
import numpy as np


def _main():
    start = time.time()

    args = get_args()

    x = load(args.x)
    y = load(args.y)

    delta_phi = get_polcal(args)

    if not (args.ds or args.t):
        print("Not generating dynamic spectra or time series - exiting!")
        sys.exit(1)

    if not (args.X or args.Y or args.I or args.Q or args.U or args.V):
        print("Not saving any polarisation or Stokes data - exiting!")
        sys.exit(2)

    if args.t:
        if args.I or args.Q or args.U or args.V:
            print("Calculating Stokes parameters")
            calculate_stokes(args, x, y, args.o, "t")

    if args.ds:
        print("Generating x and y dynamic spectra")
        x_ds = save_load("TEMP_x_ds.npy", generate_dynspec(x))
        y_ds = save_load("TEMP_y_ds.npy", generate_dynspec(y))

        if args.X:
            save(x_ds, args.o, "x", "dynspec")

        if args.Y:
            save(y_ds, args.o, "y", "dynspec")

        if args.I or args.Q or args.U or args.V:
            print("Calculating Stoke parameters")
            stokes_fnames = calculate_stokes(
                args, x_ds, y_ds, args.o, "dynspec", delta_phi=delta_phi
            )
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
    parser.add_argument(
        "-p", 
        type=str,
        default=None,
        help="Polarisation calibratrion solutions",
    )
    parser.add_argument(
        "-f", type=float, help="Central frequency (MHz)"
    )
    parser.add_argument(
        "--bw", type=float, help="Bandwidth (MHz)"
    )
    return parser.parse_args()


def load(fname):
    print(f"Loading {fname}")
    return np.load(fname, mmap_mode="r")


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


def calculate_stokes(args, x, y, outfile, type, delta_phi=None):
    # lambda functions for each of the Stokes parameters
    stokes = {
        "I": lambda x, y: np.abs(x) ** 2 + np.abs(y) ** 2,
        "Q": lambda x, y: np.abs(x) ** 2 - np.abs(y) ** 2,
        "U": lambda x, y: 2 * np.real(np.conj(x) * y),
        "V": lambda x, y: 2 * np.imag(np.conj(x) * y),
    }

    stk_args = [args.I, args.Q, args.U, args.V]
    fnames = []

    #stks = ["I"] if type == "t" else ["I", "Q", "U", "V"]
    stks = ["I", "Q", "U", "V"]
    pars = []

    for idx, stk in enumerate(stks):
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
            
            pars.append(save_load(f"TEMP_{stk}_{type}.npy", par))

    if delta_phi is not None:
        print("Applying polarisation calibration solutions")
        delta_phi_rpt = np.repeat(delta_phi[:, np.newaxis], x.shape[0], axis=1)        
        delta_phi_rpt = save_load("TEMP_delta_phi_rpt.npy", delta_phi_rpt)

        cos_phi = save_load("TEMP_cos_phi.npy", np.cos(delta_phi_rpt))
        sin_phi = save_load("TEMP_sin_phi.npy", np.sin(delta_phi_rpt))

        U_prime = pars[2]
        V_prime = pars[3]

        # apply polcal solutions via rotation matrix
        print("U")
        U = U_prime * cos_phi - V_prime * sin_phi
        print("V")
        V = U_prime * sin_phi + V_prime * cos_phi

        pars = [pars[0], pars[1], U, V]

    for i, par in enumerate(pars):
        fnames.append(save(par, outfile, stks[i], type))
    
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


def get_polcal(args):
    """Create Delta phi(f) array for polarisation calibration of Stokes
    U and V dynamic spectra.

    :param args: Command line arguments
    :type args: Namespace
    :return delta_phi: Frequency-dependent polcal solution
    :type delta_phi: np.ndarray
    """
    if args.p is not None:
        freqs = np.arange(
            args.f - args.bw/2,
            args.f + args.bw/2
        ) + 0.5
        
        lines = [float(line.rstrip("\n")) for line in open(args.p)]
        # lines == [delay, offset]
        ply = np.poly1d(lines)
        delta_phi = ply(freqs)
    else:
        delta_phi = None

    return delta_phi


def save_load(fname, x):
    """Save an array and reload it as a memory map to reduce memory
    usage

    :param fname: Filename to save to
    :type fname: str
    :param x: Array to save and load as memory map
    :type x: np.ndarray
    :return: Memory map of x
    :rtype: np.memmap
    """
    print(f"Saveloading {fname}")
    np.save(fname, x)
    del x
    return np.load(fname, mmap_mode="r")


if __name__ == "__main__":
    _main()
