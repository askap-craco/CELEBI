import sys
import time

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
        
        if args.p:
            plot(args, stokes_fnames)

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
        "-p", action="store_true", default=False, help="Generate plots"
    )
    parser.add_argument(
        "-f", type=float, help="Central frequency (MHz). Required for plotting."
    )
    parser.add_argument(
        "-l", "--label", type=str, help="FRB label. For saving plots."
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
        "i": lambda x, y: np.abs(x) ** 2 + np.abs(y) ** 2,
        "q": lambda x, y: np.abs(x) ** 2 - np.abs(y) ** 2,
        "u": lambda x, y: 2 * np.real(np.conj(x) * y),
        "v": lambda x, y: 2 * np.imag(np.conj(x) * y),
    }

    stk_args = [args.I, args.Q, args.U, args.V]
    fnames = []

    for idx, stk in enumerate(["i", "q", "u", "v"]):
        if stk_args[idx]:
            print(f"Calculating {stk} {type}")
            par = stokes[stk](x, y)
            if type == "dynspec":
                par_norm = normalise(par)
                del par
                par = par_norm.transpose()
                del par_norm

            fnames.append(save(par, outfile, stk, type))
    
    return fnames


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


def reduce(a, n, axis=0):
    """Reduces the time resolution of a given array by a factor n.

    Parameters
    ----------
    a : array
        input array to be reduced

    n : int
        factor to reduce by

    axis : int
        axis along which

    Returns
    -------
    array
        reduced array
    """

    if n > 1:
        if axis == 1:
            a = a.transpose()
        a_red = []
        for i in range(int(a.shape[0] / n)):
            a_red.append(np.sum(a[i * n : (i + 1) * n], axis=0))

        a_red = np.array(a_red) / np.sqrt(n)

        if axis == 1:
            a_red = a_red.transpose()

        return np.array(a_red)
    else:
        return a


def plot_IQUV_dts(
    ds,
    f0,
    facs=[10, 30, 100, 300],
    time_range=100,
    peaks=None,
    fig=None,
    axs=None,
    title=True,
    xlabel=True,
):
    new_peaks = []

    if axs is None or fig is None:
        fig, axs = plt.subplots(nrows=1, ncols=len(facs), figsize=(12, 6))

    for i, dt in enumerate(facs):
        ds_red = reduce(ds, dt, axis=1)
        peak = peaks[i] if peaks else np.argmax(np.sum(ds_red, axis=0))
        new_peaks.append(peak)
        plot_ds = ds_red[:, peak - time_range : peak + time_range]
        ax = axs.flatten()[i]
        ax.imshow(
            plot_ds,
            aspect="auto",
            interpolation="none",
            extent=(
                -time_range * dt / 1e3,
                time_range * dt / 1e3,
                f0 - 336 / 2,
                f0 + 336 / 2,
            ),
        )

        if xlabel:
            ax.set_xlabel("Time (ms)")

        if title:
            ax.set_title(f"dt = {dt} us")

        if i == 0:
            ax.set_ylabel("Frequency (MHz)")
        else:
            ax.set_yticks([])

    return (fig, axs, new_peaks)


def plot(args, stokes_fnames):
    stks = [np.load(f, mmap_mode="r") for f in stokes_fnames]

    # IQUV over four timescales
    fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(10, 10))
    for i, ds in enumerate(stks):
        fig, a, peaks = plot_IQUV_dts(
            ds,
            args.f,
            time_range=50,
            fig=fig,
            axs=axs[i],
            xlabel=True if i == 3 else False,
            title=True if i == 0 else False,
            peaks=None if i == 0 else peaks,
        )
    plt.tight_layout()
    plt.savefig(f"{args.label}_IQUV_dts.png")

if __name__ == "__main__":
    _main()
