import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np


def _main():
    args = get_args()

    with open(args.s) as fnames_file:
        fnames = fnames_file.readlines()
        fnames = [fname.strip() for fname in fnames]
    
    plot(args, fnames)


def get_args():
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

    parser = ArgumentParser(
        description="Plots FRB Stokes dynamic spectra"
    )
    parser.add_argument(
        "-s", type=str, required=True, help="File containing names of dynspec files"
    )
    parser.add_argument(
        "-f", type=float, help="Central frequency (MHz)"
    )
    parser.add_argument(
        "-l", "--label", type=str, help="FRB label"
    )
    parser.add_argument(
        "-d", "--DM", type=str, help="Dedispersion DM"
    )

    return parser.parse_args()


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
    facs=[1, 3, 10, 30, 100, 300],
    time_range=100,
    labels=None,
    fname=None,
):
    fig = plt.figure(figsize=(2.5 * len(facs), 2.5 * len(ds_list)))
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
    plot_IQUV_dts(stks, args.f, labels=["I", "Q", "U", "V"], fname=f"{args.label}_IQUV_dts_{args.DM}.png")


if __name__ == "__main__":
    _main()