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
    ds_list,
    f0,
    facs=[10, 30, 100, 300],
    time_range=100,
    labels=None,
    fname=None,
):
    fig = plt.figure(figsize=(2.5 * len(ds_list), 2.5 * len(facs)))
    spec = gs.GridSpec(nrows=len(ds_list) + 1, ncols=3 * len(facs) + 1)

    peaks = []

    prof_axs = [fig.add_subplot(spec[0, 3 * j : 3 * (j + 1)]) for j in range(len(facs))]

    for i, ds in enumerate(ds_list):
        for j, dt in enumerate(facs):
            ds_red = reduce(ds, dt, axis=1)
            if len(peaks) == j:
                # first time through timescales
                peaks.append(np.argmax(np.sum(ds_red, axis=0)))

            peak = peaks[j]

            plot_ds = ds_red[:, peak - time_range : peak + time_range]
            plot_ts = np.sum(plot_ds, axis=0)
            t = np.linspace(
                -time_range * dt / 1e3, time_range * dt / 1e3, plot_ts.shape[0]
            )

            prof_ax = prof_axs[j]
            ds_ax = fig.add_subplot(spec[i + 1, 3 * j : 3 * (j + 1)])

            prof_ax.step(t, plot_ts, lw=1, label=labels[i])
            prof_ax.set_xlim(t[0], t[-1])
            ds_ax.imshow(
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

            prof_ax.set_xticks([])
            if i == 0:
                # first row
                prof_ax.set_title(f"dt = {dt} us")
                ds_ax.set_xticks([])
            elif i == len(ds_list) - 1:
                # last row
                ds_ax.set_xlabel("Time (ms)")
            else:
                ds_ax.set_xticks([])

            prof_ax.set_yticks([])
            if j == 0:
                # first column
                ds_ax.set_ylabel("Frequency (MHz)")
            else:
                ds_ax.set_yticks([])

        # spectrum column
        plot_spec = ds_red[:, peak]
        freqs = np.linspace(f0 - 336 / 2, f0 + 336 / 2, plot_spec.shape[0])

        spec_ax = fig.add_subplot(spec[i + 1, -1])
        spec_ax.step(plot_spec, freqs[::-1], c="k", lw=1)
        spec_ax.axvline(0, c="k", lw=1, ls="--")

        spec_ax.set_yticks([])
        spec_ax.set_xticks([])

        spec_ax.set_ylim(freqs[0], freqs[-1])

        spec_ax.yaxis.set_label_position("right")
        spec_ax.set_ylabel(labels[i], rotation=0, fontsize=20, labelpad=20)

    prof_ax.legend(bbox_to_anchor=(1, 0.5), loc="center left")

    plt.subplots_adjust(
        wspace=0, hspace=0, left=0.075, right=0.95, top=0.975, bottom=0.05
    )

    if fname:
        plt.savefig(fname)


def plot(args, stokes_fnames):
    stks = [np.load(f, mmap_mode="r") for f in stokes_fnames]

    # IQUV over four timescales
    plot_IQUV_dts(stks, args.f0, labels=["I", "Q", "U", "V"], fname=f"{args.label}_IQUV_dts.png")

if __name__ == "__main__":
    _main()
