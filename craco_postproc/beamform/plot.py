import matplotlib as mpl
mpl.use('agg')

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

    abs_peak = None
    peaks = []

    prof_axs = [fig.add_subplot(spec[0, 3 * j : 3 * (j + 1)]) for j in range(len(facs))]

    for i, ds in enumerate(ds_list):
        print(f"{i}")
        j = len(prof_axs) - 1
        prev_fac = None
        for k, dt in enumerate(facs[::-1]):
            print(f"  {j}, {dt}")
            ds_red = reduce(ds, dt, axis=1)
            ts_red = np.sum(ds_red, axis=0)
            if abs_peak is None:
                abs_peak = (np.argmax(ts_red) * dt,
                            np.argmax(ts_red) * dt+dt+1)
            
            if len(peaks) < len(prof_axs):
                peaks.append(
                    abs_peak[0] // dt + 
                    np.argmax(ts_red[abs_peak[0]//dt:abs_peak[1]//dt])
                )

            # choose peak within abs_peak
            peak = peaks[k]

            print(f"  peak = {peak}")
            
            plot_ds = ds_red[
                :, 
                max(0, peak - time_range) : min(ds_red.shape[1], peak + time_range)
            ]

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

            j -= 1
            prev_fac = dt

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

    return peaks[-1]  # 1 us peak


def crop(stks, peak, facs, labels, time_range, prefix):
    c = slice(peak - time_range//2, peak + time_range//2, 1)
    for i, dt in enumerate(facs):
        for j, s in enumerate(stks):
            s_red = reduce(s, dt, axis=1)
            red_c = slice(
                c.start // dt,
                c.stop // dt,
                1
            )
            np.save(f"{prefix}{dt:03d}us_{labels[j]}.npy", s_red[:, red_c])


def plot_ts(ds, fname):
    ts = np.sum(ds, axis=0)
    ts = reduce(ts, 1000)
    plt.figure(figsize=(10, 5))
    plt.plot(ts)
    plt.xlabel("Time (ms)")
    plt.ylabel("I")
    plt.tight_layout()
    plt.savefig(fname)


def plot(args, stokes_fnames):
    stks = [np.load(f, mmap_mode="r") for f in stokes_fnames]

    # IQUV over four timescales
    peak = plot_IQUV_dts(
        stks, 
        args.f,
        facs=[1, 3, 10, 30, 100, 300],
        labels=["I", "Q", "U", "V"], 
        fname=f"{args.label}_IQUV_dts_{args.DM}.png"
    )

    crop(
        stks,
        peak,
        [1, 3, 10, 30, 100, 300],
        ["I", "Q", "U", "V"],
        100*1000,    # 100 ms
        f"crops/{args.label}_{args.DM}_"
    )

    # plot full Stokes I time series at 1 ms time resolution
    plot_ts(stks[0], f"{args.label}_I_{args.DM}.png")


if __name__ == "__main__":
    _main()
