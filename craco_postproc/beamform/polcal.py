import argparse
from functools import reduce

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as un
from scipy.optimize import bisect, curve_fit


def _main():
    args = get_args()

    print("Loading IQUV")
    I = np.load(args.i, mmap_mode="r")
    Q = np.load(args.q, mmap_mode="r")
    U = np.load(args.u, mmap_mode="r")
    V = np.load(args.v, mmap_mode="r")

    freqs, df, dt = get_freqs(
        args.f * un.MHz, args.b * un.MHz, I.shape[0], args.reduce_df
    )

    period = args.p * un.s
    print(f"period = {period}")

    # get Stokes spectra and noise
    print("Getting spectra")
    I_f, I_noise, peak = get_pulsar_spec(I, period, dt, args, "I")
    Q_f, Q_noise, peak = get_pulsar_spec(Q, period, dt, args, "Q", peak=peak)
    U_f, U_noise, peak = get_pulsar_spec(U, period, dt, args, "U", peak=peak)
    V_f, V_noise, peak = get_pulsar_spec(V, period, dt, args, "V", peak=peak)

    I_noisesub = I_f - I_noise
    Q_noisesub = Q_f - Q_noise
    U_noisesub = U_f - U_noise
    V_noisesub = V_f - V_noise

    S_noisesub = [I_noisesub, Q_noisesub, U_noisesub, V_noisesub]
    S_names = ["I", "Q", "U", "V"]

    if args.plot:
        # plot Stokes
        print("Plotting Stokes")
        fig, axs = plt.subplots(2, 2)
        axs = axs.flatten()
        for i in range(4):
            axs[i] = ax_plot(axs[i], freqs, S_noisesub[i])
            axs[i].set_title(f"Stokes {S_names[i]}")
            if i < 2:
                axs[i].set_xticks([])
            else:
                axs[i].set_xlabel("Frequency (MHz)")

            if i % 2 == 0:
                axs[i].set_ylabel("Intensity (arb. units)")
        plt.tight_layout()
        fig.savefig(f"{args.plotdir}/{args.label}_Stokes.png")
        plt.close(fig)

    S_ratio = [
        Q_noisesub / I_noisesub,
        U_noisesub / I_noisesub,
        V_noisesub / I_noisesub,
    ]

    if args.plot:
        # plot Stokes ratios
        print("Plotting Stokes ratios")
        fig, axs = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)
        for i in range(3):
            axs[i] = ax_plot(
                axs[i], freqs, S_ratio[i], xlabel="Frequency (MHz)"
            )
            axs[i].set_title(f"{S_names[i+1]}/{S_names[0]}")
        axs[0].set_ylim(-1, 1)
        plt.tight_layout()
        fig.savefig(f"{args.plotdir}/{args.label}_Stokes_ratios.png")
        plt.close(fig)

    print("Calculating and fitting polarisation fraction")
    pol_f = S_ratio[0] ** 2 + S_ratio[1] ** 2 + S_ratio[2] ** 2
    pol_fit = np.polyfit(freqs, pol_f, 0)[0]
    print(f"Pol frac fit: {pol_fit}")

    if args.plot:
        fig, ax = plt.subplots()
        ax = ax_plot(
            ax,
            freqs,
            pol_f,
            xlabel="Frequency (MHz)",
            ylabel="Polarisation fraction",
        )
        ax.axhline(pol_fit)
        plt.tight_layout()
        fig.savefig(f"{args.plotdir}/{args.label}_polfrac.png")

    print("Calculating Polarisation Angle")
    pol_ang = np.arctan2(U_noisesub, Q_noisesub) / 2

    if args.plot:
        fig, ax = plt.subplots()
        ax = ax_plot(
            ax,
            freqs,
            pol_ang,
            type="scatter",
            xlabel="Frequency (MHz)",
            ylabel="Pol angle",
        )
        ax.set_ylim(-np.pi / 2, np.pi / 2)
        plt.yticks(
            np.linspace(-1 / 2, 1 / 2, 5) * np.pi,
            [
                r"-$\frac{\pi}{2}$",
                r"-$\frac{\pi}{4}$",
                r"$0$",
                r"$\frac{\pi}{4}$",
                r"$\frac{\pi}{2}$",
            ],
        )
        plt.tight_layout()
        fig.savefig(f"{args.plotdir}/{args.label}_polang.png")

    # rm, offset, stokes_ratio_pks = polcal_ref(args)
    rm = 37.43534317649634
    offset = -1.80339213320941

    # determine ASKAP psi_sky
    # psi'(nu) = psi(nu) + psi_sky
    # TODO: tidy this up, copied almost verbatim from polcal_with_VELA_for_shivani.ipynb
    def QoverI_askap(pa, psi_sky, Lamp):
        return Lamp * np.cos(2 * (pa + psi_sky))

    pa_askap = faraday_angle(freqs.value, rm, offset)
    Q_askap = np.copy(S_noisesub[1] / S_noisesub[0])
    popt, pcov = curve_fit(QoverI_askap, pa_askap, Q_askap, p0=[-0.8, 0.95])
    psi_sky = popt[0]
    L_amp = popt[1]

    if args.plot:
        fig, ax = plt.subplots()
        ax = ax_plot(ax, freqs, Q_askap, label=r"$Q/I$ (askap)$")
        ax = ax_plot(
            ax,
            freqs,
            QoverI_askap(pa_askap, *popt),
            label=r"$\psi_{sky}$ fit",
            c="r",
            type="line",
        )
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Q/I")
        ax.legend()
        plt.tight_layout()
        fig.savefig(f"{args.plotdir}/{args.label}_QI_fit.png")

    # Solve using matrices
    tan_phi = np.full((freqs.shape[0], 2), np.nan)
    u_prime = S_noisesub[2] / S_noisesub[0]  # U/I
    v_prime = S_noisesub[3] / S_noisesub[0]  # V/I
    for i, f in enumerate(freqs.value):
        pa_prime = psi_sky + faraday_angle(f, rm, offset)

        u = L_amp * np.sin(2 * pa_prime)
        v = -0.05

        A = np.array([[u + v, v - u],
                      [u - v, u + v]])
        x = np.array([[u_prime[i] + v_prime[i]],
                      [u_prime[i] - v_prime[i]]])
        y = np.matmul(np.linalg.inv(A), x)

        tan_phi[i, 0] = y[0]
        tan_phi[i, 1] = y[1]

    phi = np.arctan2(tan_phi[:, 1], tan_phi[:, 0])
    pfit, pcov = np.polyfit(freqs.value, phi, 1, cov=True)
    perr = np.sqrt(np.diag(pcov))

    delay_ns = pfit[0] / (2 * np.pi * 1e6)
    delay_err_ns = perr[0] / (2 * np.pi * 1e6)
    offset = pfit[1]
    offset_err = perr[1]

    print(
        f"Delay:  {delay_ns} +- {delay_err_ns} s "
        f"({int(delay_err_ns/delay_ns*100)}% error)"
    )
    print(
        f"Offset: {offset} +- {offset_err} "
        f"({int(offset_err/offset*100)}% error)"
    )

    if args.plot:
        fig, ax = plt.subplots()
        ax = ax_plot(ax, freqs, phi, label="Solved")
        p = np.poly1d([delay_ns * (2 * np.pi * 1e6), offset])
        ax = ax_plot(
            ax,
            freqs,
            p(freqs.value),
            label="Fit",
            c="r",
            type="line",
        )
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel(r"$\Delta \Phi$")
        ax.legend()
        plt.tight_layout()
        fig.savefig(f"{args.plotdir}/{args.label}_solved_deltaphi.png")


    with open(args.o, "w") as f:
        f.write(f"delay\t{delay_ns}\t+-\t{delay_err_ns}\n")
        f.write(f"offset\t{offset}\t+-\t{offset_err}")


def get_args() -> argparse.Namespace:
    """Parse command line arguments

    :return: Command line argument paramters
    :rtype: :class:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser(
        description="Determine and apply polarisation "
        "calibration solutions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", type=str, help="Stokes I dynamic spectrum")
    parser.add_argument("-q", type=str, help="Stokes Q dynamic spectrum")
    parser.add_argument("-u", type=str, help="Stokes U dynamic spectrum")
    parser.add_argument("-v", type=str, help="Stokes V dynamic spectrum")
    parser.add_argument("-p", type=float, help="Pulsar period in s")
    parser.add_argument("-f", type=float, help="Centre frequency in MHz")
    parser.add_argument("-b", type=float, help="Bandwidth in MHz")
    parser.add_argument("-l", "--label", type=str, help="Label for plots")
    parser.add_argument("-r", "--ref", type=str, help="Reference Parkes data")
    parser.add_argument("-o", type=str, help="File to output delay and offset")
    parser.add_argument(
        "--reduce_df",
        type=int,
        default=1,
        help="If given, reduces the frequency resolution by the given factor.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        default="False",
        help="If given, will plot at various stages",
    )
    parser.add_argument("--plotdir", type=str, help="Directory to plot into.")
    parser.add_argument(
        "--refplotdir",
        type=str,
        help="Directory to plot reference data into.",
    )
    return parser.parse_args()


def get_freqs(c_freq, bw, n_chan, reduce_df):
    # "tuple[np.ndarray, astropy.Quantity, astropy.Quantity]"
    """Determine fine channel frequencies, and frequency & time
    resolutions.

    :param c_freq: Central frequency
    :type c_freq: :class:`astropy.Quantity`
    :param bw: Bandwidth
    :type bw: :class:`astropy.Quantity`
    :param n_chan: Number of fine channels
    :type n_chan: int
    :param reduce_df: Factor to reduce frequency resolution by
    :type reduce_df: int
    :return: Frequency array, frequency resolution, time resolution
    :rtype: tuple[np.ndarray, astropy.Quantity, astropy.Quantity]
    """
    print("Determining freqs")
    freqs = np.linspace(
        c_freq - bw / 2, c_freq + bw / 2, n_chan, endpoint=False
    )
    freqs = freqs[:: reduce_df]

    df = freqs[1] - freqs[0]
    dt = (1 / df) * reduce_df
    print(f"df = {df}")
    print(f"dt = {dt.to(un.us)}")

    return freqs, df, dt


def fold_ds(a, p, dt):
    """
    Fold a along last axis according to period p and time resolution dt
    """
    t_max = a.shape[-1] * dt.to(un.s)
    dt_per_p = int(p / dt)
    num_periods = int(t_max / p)

    print(f"{num_periods} periods")

    new_t_ax_shape = dt_per_p
    fold_arr = np.zeros((a.shape[0], new_t_ax_shape, num_periods))

    for i in range(num_periods):
        fold_arr[:, :, i] = a[:, i * dt_per_p : (i + 1) * dt_per_p]

    return np.sum(fold_arr, axis=2)


def reduce(a, n, transpose=False, pbar=True):
    """Reduces the time resolution of a given array by a factor n.

    Reduces along the 0th axis, make sure this one is the time axis!
    Returned array is the time-averaged version of the imported array.

    Parameters
    ----------
    a : array
        input array to be reduced

    n : int
        factor to reduce by

    transpose : bool
        does the function need to transpose the array before reducing?
        Set to True iff axis 0 is not time
        (Default value = False)

    pbar : bool, optional

    Returns
    -------
    array
        reduced array
    """

    if n > 1:
        a = a.transpose() if transpose else a
        A_red = []
        for i in range(int(a.shape[0] / n)):
            A_red.append(np.sum(a[i * n : (i + 1) * n], axis=0))
        return (
            np.array(A_red) / np.sqrt(n)
            if transpose
            else np.array(A_red) / np.sqrt(n)
        )
    else:
        return a


def get_pulsar_spec(
    ds, period, dt, args, par, peak=None, red_fac=100, width=10
):
    folded_ds = fold_ds(ds, period, dt)
    folded_ds_red = reduce(folded_ds, red_fac, transpose=True).transpose()

    if args.reduce_df > 1:
        folded_ds_red = reduce(folded_ds_red, args.reduce_df)

    peak = peak if peak else np.argmax(np.sum(folded_ds_red, axis=0))

    turn_on = peak - width
    turn_off = peak + width

    fmin = args.f - args.b / 2
    fmax = args.f + args.b / 2
    if args.plot:
        plt.imshow(
            folded_ds_red,
            extent=(0, folded_ds_red.shape[1], fmin, fmax),
            aspect="auto",
            interpolation="none",
        )
        plt.axvline(peak, c="r", lw=1)
        plt.axvline(turn_on, c="r", linestyle="--", lw=1)
        plt.axvline(turn_off, c="r", linestyle="--", lw=1)
        plt.xlabel("Time")
        plt.ylabel("Frequency (MHz)")
        plt.title(f"Stokes {par} (on peak)")
        plt.tight_layout()
        plt.savefig(f"{args.plotdir}/{args.label}_Stokes{par}_peak.png")

    spec = np.sum(folded_ds_red[:, turn_on:turn_off], axis=1)

    # roll folded data halfway along time axis and get noise with same indexes
    folded_ds_red = np.roll(
        folded_ds_red, int(folded_ds_red.shape[1] / 2), axis=1
    )

    if args.plot:
        plt.imshow(
            folded_ds_red,
            extent=(0, folded_ds_red.shape[1], fmin, fmax),
            aspect="auto",
            interpolation="none",
        )
        plt.axvline(peak, c="r", lw=1)
        plt.axvline(turn_on, c="r", linestyle="--", lw=1)
        plt.axvline(turn_off, c="r", linestyle="--", lw=1)
        plt.xlabel("Time")
        plt.ylabel("Frequency (MHz)")
        plt.title(f"Stokes {par} (off peak)")
        plt.tight_layout()
        plt.savefig(f"{args.plotdir}/{args.label}_Stokes{par}_noise.png")

    noise = np.sum(folded_ds_red[:, turn_on:turn_off], axis=1)

    return spec, noise, peak


# Parkes PA fit psi(nu): Get RM and offset
def faraday_angle(freq_mhz, RM=30, offset=0):
    lamb = 3e8 / (freq_mhz * 1e6)
    return RM * np.power(lamb, 2) + offset


def polcal_ref(args):
    plotdir = args.refplotdir
    vela_pks = np.loadtxt(
        args.ref,
        dtype={
            "names": ("freq_ghz", "I", "Q", "U", "V", "PA"),
            "formats": ("f4", "f4", "f4", "f4", "f4", "f4"),
        },
    )
    freq_mhz = (vela_pks["freq_ghz"] * un.GHz).to(un.MHz)
    stokes_ratio_pks = [
        np.copy(vela_pks["Q"]) / vela_pks["I"],
        np.copy(vela_pks["U"]) / vela_pks["I"],
        np.copy(vela_pks["V"]) / vela_pks["I"],
    ]
    stokes_ratio_name = ["Q/I", "U/I", "V/I"]
    pa = np.copy(vela_pks["PA"])

    if args.plot:
        fig, ax = plt.subplots()
        cs = ["r", "g", "b"]
        for i in range(3):
            ax = ax_plot(
                ax,
                freq_mhz,
                stokes_ratio_pks[i],
                c=cs[i],
                label=stokes_ratio_name[i],
                xlabel="Frequency (MHz)",
            )
        ax.legend()
        plt.tight_layout()
        fig.savefig(f"{plotdir}/vela_parkes_Stokes_ratios.png")

    if args.plot:
        fig, ax = plt.subplots()
        ax = ax_plot(
            ax, freq_mhz, pa, xlabel="Frequency (MHz)", ylabel="Pol angle"
        )
        plt.tight_layout()
        fig.savefig(f"{plotdir}/vela_parkes_polang.png")

    # determine p for Parkes
    linearIntensity = (
        stokes_ratio_pks[0] ** 2
        + stokes_ratio_pks[1] ** 2
        + stokes_ratio_pks[2] ** 2
    )
    nan = np.isnan(linearIntensity)
    np.mean(linearIntensity[~nan])
    if args.plot:
        fig, ax = plt.subplots()
        ax = ax_plot(
            ax,
            freq_mhz,
            linearIntensity,
            xlabel="Frequency (MHz)",
            ylabel="Polarisation fraction",
        )
        plt.tight_layout()
        fig.savefig(f"{plotdir}/vela_parkes_polfrac.png")

    nonzero = np.where(pa != 0)
    pa_tofit = np.copy(pa)
    popt, pcov = curve_fit(
        faraday_angle,
        freq_mhz[nonzero].value,
        pa_tofit[nonzero] / 180.0 * np.pi,
        p0=(-30, -0.2),
    )
    rm = popt[0]
    offset = popt[1]
    print(("RM", rm, "offset", offset))
    print(pcov)

    if args.plot:
        fig, ax = plt.subplots()
        ax = ax_plot(
            ax,
            freq_mhz[nonzero],
            pa_tofit[nonzero] / 180.0 * np.pi,
            label="Parkes",
        )
        ax = ax_plot(
            ax,
            freq_mhz,
            faraday_angle(freq_mhz.value, popt[0], popt[1]),
            label="PA fit",
            c="r",
            type="line",
        )
        ax.set_ylabel("(rad)")
        ax.set_xlabel("Frequency (MHz)")
        ax.legend()
        plt.tight_layout()
        fig.savefig(f"{plotdir}/vela_parkes_polang_fit.png")

    # Determine Parkes L/I (L_amp) via fitting for Q/I using psi(nu) from above
    pa = faraday_angle(freq_mhz.value, rm, offset)

    nan = np.isnan(stokes_ratio_pks[0])

    def QoverI(pa, amp):
        return amp * np.cos(2 * pa)

    def UoverI(pa, amp):
        return amp * np.sin(2 * pa)

    return rm, offset, stokes_ratio_pks


def ax_plot(
    ax,
    x,
    y,
    xlabel=None,
    ylabel=None,
    label=None,
    c="k",
    lw=1,
    type="step",
    **kwargs,
):
    """
    Plot y vs x on ax.

    Type can be one of:
        step
        line
        scatter
    """
    if type == "step":
        ax.step(x, y, where="mid", c=c, lw=lw, label=label, **kwargs)
    elif type == "line":
        ax.plot(x, y, c=c, lw=lw, label=label, **kwargs)
    elif type == "scatter":
        ax.scatter(x, y, c=c, label=label, s=5, **kwargs)
    else:
        print(f"Plot type {type} not valid")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax


if __name__ == "__main__":
    _main()
