import argparse
from functools import reduce

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as un
from scipy.optimize import curve_fit


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

    t = np.arange(I_f.shape[0])

    I_noisesub = I_f - I_noise
    Q_noisesub = Q_f - Q_noise
    U_noisesub = U_f - U_noise
    V_noisesub = V_f - V_noise

    S_noisesub = [I_noisesub, Q_noisesub, U_noisesub, V_noisesub]
    S_names = ["I", "Q", "U", "V"]

    S_ratio = np.array([
        Q_noisesub / I_noisesub,
        U_noisesub / I_noisesub,
        V_noisesub / I_noisesub,
    ])

    if args.nopeak:  # fit in each time step
        pol_fs = np.zeros((t.shape[0], freqs.shape[0]))
        pol_fits = np.zeros(t.shape)
        pas = np.zeros(pol_fs.shape)
        pa_fits = np.zeros(pol_fs.shape)
        Q_askaps = np.zeros(pol_fs.shape)
        L_amps = np.zeros(pol_fits.shape)
        psi_skys = np.zeros(pol_fits.shape)
        Q_fits = np.zeros(pol_fs.shape)
        phis = np.zeros(pol_fs.shape)
        delays = np.zeros(pol_fits.shape)
        offsets = np.zeros(pol_fits.shape)
        phi_fits = np.zeros(pol_fs.shape)

        for i in range(t.shape[0]):
            print(f"{i+1}/{t.shape[0]}")
            pol_fs[i], pol_fits[i] = fit_pol_frac(freqs, S_ratio[:, i])

            pas[i], pa_fits[i] = fit_pol_ang(
                freqs, U_noisesub[i], Q_noisesub[i]
            )

            Q_askaps[i], L_amp_fit, psi_sky_fit, Q_fits[i] = fit_psi_sky(
                freqs, S_ratio[:, i], pa_fits[i]
            )
            L_amps[i] = L_amp_fit[0]
            psi_skys[i] = psi_sky_fit[0]

            phis[i], delay_fit, offset_fit, phi_fits[i] = fit_phi(
                freqs, S_ratio[:, i], pa_fits[i], L_amp_fit[0], psi_sky_fit[0]
            )
            delays[i] = delay_fit[0]
            offsets[i] = offset_fit[0]
        
        if args.plot:
            ext = (
                t[0], t[-1],
                freqs[0].value, freqs[-1].value
            )
            # Stokes
            fig, axs = plt.subplots(2, 2)
            axs = axs.flatten()
            for i in range(4):
                axs[i] = im_plot(
                    axs[i], 
                    S_noisesub[i], 
                    title=f"Stokes {S_names[i]}",
                    extent=ext,
                )

                if i < 2:
                    axs[i].set_xticks([])
                else:
                    axs[i].set_xlabel("Time (samp.)")

                if i % 2 == 0:
                    axs[i].set_ylabel("Frequency (MHz)")
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_Stokes.png")

            # Stokes ratios
            fig, axs = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)
            for i in range(3):
                axs[i] = im_plot(
                    axs[i], 
                    S_ratio[i], 
                    xlabel="Time (samp.)", 
                    title=f"{S_names[i+1]}/{S_names[0]}",
                    vmin=-1, 
                    vmax=1,
                    extent=ext,
                )
            axs[0].set_ylabel("Frequency (MHz)")
            fig.savefig(f"{args.plotdir}/{args.label}_Stokes_ratios.png")

            # polarisation fraction
            fig, axs = plt.subplots(
                nrows=2, ncols=1, figsize=(8, 12), sharex=True
            )
            ax = ax_plot(
                axs[0],
                t,
                pol_fits,
                ylabel="Fit polarisation fraction"
            )
            axs[0].set_ylim(0, 1)
            ax = im_plot(
                axs[1],
                pol_fs,
                xlabel="Time (samp.)",
                ylabel="Frequency (MHz)",
                extent=ext,
                vmax=1,
                cbar="Polarisation fraction",
                orientation="horizontal",
            )
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_polfrac.png")

            # polarisation angle
            fig, axs = plt.subplots(
                nrows=2, ncols=1, figsize=(8,12), sharex=True, sharey=True
            )
            axs[0] = im_plot(
                axs[0],
                pa_fits,
                ylabel="Frequency (MHz)",
                title="Fit PA",
                extent=ext,
                vmin=-np.pi/2,
                vmax=np.pi/2,
            )
            axs[1] = im_plot(
                axs[1],
                pas,
                xlabel="Time (samp.)",
                ylabel="Frequency (MHz)",
                title="Measured PA",
                extent=ext,
                vmin=-np.pi/2,
                vmax=np.pi/2,
                #cbar="PA",
                #orientation="horizontal",
            )
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_polang.png")

            # psi_sky fit
            pa_ext = (
                t[0], t[-1],
                np.min(pa_fits), np.max(pa_fits)
            )
            fig, axs = plt.subplots(
                nrows=2, ncols=1, figsize=(8,12), sharex=True, sharey=True
            )
            axs[0] = im_plot(
                axs[0],
                Q_fits,
                ylabel="PA (rad.)",
                title="Fit Q/I",
                extent=pa_ext,
                vmin=-1,
                vmax=1,
            )
            axs[1] = im_plot(
                axs[1],
                Q_askaps,
                xlabel="Time (samp.)",
                ylabel="PA (rad.)",
                title="Measured Q/I",
                extent=pa_ext,
                vmin=-1,
                vmax=1,
            )
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_QI_fit.png")

            # delta phi fit
            fig, axs = plt.subplots(
                nrows=2, ncols=1, figsize=(8,12), sharex=True, sharey=True
            )
            axs[0] = im_plot(
                axs[0],
                phi_fits,
                ylabel="Frequency (MHz)",
                title=r"Fit $\Delta\phi$",
                extent=ext,
            )
            axs[1] = im_plot(
                axs[1],
                phis,
                xlabel="Time (samp.)",
                ylabel="Frequency (MHz)",
                title=r"Measured $\Delta\phi$",
                extent=ext,
            )
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_solved_deltaphi.png")

    else:
        pol_f, pol_fit = fit_pol_frac(freqs, S_ratio)

        pa, pa_fit = fit_pol_ang(freqs, U_noisesub, Q_noisesub)

        Q_askap, L_amp_fit, psi_sky_fit, Q_fit = fit_psi_sky(freqs, S_ratio, pa_fit)

        phi, delay_fit, offset_fit, phi_fit = fit_phi(
            freqs, S_ratio, pa_fit, L_amp_fit[0], psi_sky_fit[0]
        )

        if args.plot:
            # Stokes
            fig, axs = plt.subplots(2, 2)
            axs = axs.flatten()
            for i in range(4):
                axs[i] = ax_plot(
                    axs[i], freqs, S_noisesub[i], title=f"Stokes {S_names[i]}"
                )

                if i < 2:
                    axs[i].set_xticks([])
                else:
                    axs[i].set_xlabel("Frequency (MHz)")

                if i % 2 == 0:
                    axs[i].set_ylabel("Intensity (arb.)")
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_Stokes.png")

            # Stokes ratios
            fig, axs = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)
            for i in range(3):
                axs[i] = ax_plot(
                    axs[i], 
                    freqs,
                    S_ratio[i], 
                    xlabel="Frequency (MHz)",
                    title=f"{S_names[i+1]}/{S_names[0]}"
                )

            axs[0].set_ylim(-1, 1)
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_Stokes_ratios.png")

            # polarisation fraction
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

            # polarisation angle
            fig, ax = plt.subplots()
            ax = ax_plot(
                ax,
                freqs,
                pa,
                type="scatter",
                xlabel="Frequency (MHz)",
                ylabel="Pol angle",
                label="PA"
            )
            ax = ax_plot(
                ax,
                freqs,
                pa_fit,
                type="line",
                xlabel="Frequency (MHz)",
                ylabel="Pol angle",
                c="r",
                label="Fit"
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

            # psi_sky fit
            fig, ax = plt.subplots()
            ax = ax_plot(
                ax,
                pa_fit,
                Q_askap,
                label=r"$Q/I$ (askap)$",
                type="scatter"
            )
            ax = ax_plot(
                ax,
                pa_fit,
                QoverI_askap(pa_fit, L_amp_fit[0], psi_sky_fit[0]),
                label=r"$\psi_{sky}$ fit",
                c="r",
                type="line",
            )
            ax.set_xlabel("PA (rad)")
            ax.set_ylabel("Q/I")
            ax.legend()
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_QI_fit.png")

            # delta phi fit
            fig, ax = plt.subplots()
            ax = ax_plot(ax, freqs, phi, label="Solved")
            # for i, s in enumerate(slices):
            # p = np.poly1d([delay_fit[0], offset_fit[0]])
            ax = ax_plot(
                ax,
                freqs,
                phi_fit,
                label="Best fit",
                c="r",
                lw = 1,
                type="line",
            )
            ax.set_xlabel("Frequency (MHz)")
            ax.set_ylabel(r"$\Delta \Phi$")
            plt.tight_layout()
            fig.savefig(f"{args.plotdir}/{args.label}_solved_deltaphi.png")


    with open(args.o, "w") as f:
        f.write(f"{delay_fit[0]}\n")
        f.write(f"{offset_fit[0]}\n")


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
        default=False,
        help="If given, will plot at various stages",
    )
    parser.add_argument("--plotdir", type=str, help="Directory to plot into.")
    parser.add_argument(
        "--refplotdir",
        type=str,
        help="Directory to plot reference data into.",
    )
    parser.add_argument(
        "--nopeak",
        action="store_true",
        default=False,
        help="DEBUG: Do not take peak spectra and fit per time step"
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

    return freqs[:-1], df, dt


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
    ds, period, dt, args, par, peak=None, red_fac=100, width=5,
):
    folded_ds = fold_ds(ds, period, dt)
    folded_ds_red = reduce(folded_ds, red_fac, transpose=True).transpose()

    print(f"dt = {dt.to(un.us) * red_fac}")

    if args.reduce_df > 1:
        folded_ds_red = reduce(folded_ds_red, args.reduce_df)

    peak = peak if peak else np.argmax(np.sum(folded_ds_red, axis=0))

    turn_on = peak - width
    turn_off = peak + width + 1

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

    # [::-1] because it will be in the wrong order otherwise
    if args.nopeak:
        spec = folded_ds_red[::-1, turn_on:turn_off][:-1]
    else:
        spec = folded_ds_red[:, peak][::-1][:-1]

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

    if args.nopeak:
        noise = folded_ds_red[::-1, turn_on:turn_off][:-1]
    else:
        noise = folded_ds_red[:, peak][::-1][:-1]

    return spec.T, noise.T, peak


def fit_pol_frac(freqs, S_ratio):
    print("Calculating and fitting polarisation fraction")
    pol_f = S_ratio[0] ** 2 + S_ratio[1] ** 2 + S_ratio[2] ** 2
    pol_fit = np.polyfit(freqs, pol_f, 0)[0]
    print(f"Pol frac fit: {pol_fit}")
    return pol_f, pol_fit


def fit_pol_ang(freqs, U, Q):
    print("Calculating Polarisation Angle")
    pol_ang = np.arctan2(U, Q) / 2

    # fit rm and offset
    popt, pcov = curve_fit(faraday_angle, freqs.value, pol_ang, p0=[30, 0])
    rm, offset = popt
    rm_err, offset_err = np.diag(pcov)

    pa_fit = faraday_angle(freqs.value, rm, offset)

    print(f"rm\t= {rm:.3f}\t+- {rm_err:.5f}")
    print(f"offset\t= {offset:.3f}\t+- {offset_err:.5f}")

    return pol_ang, pa_fit


# Parkes PA fit psi(nu): Get rm and offset
def faraday_angle(freq_mhz, rm=30, offset=0):
    lamb = 3e8 / (freq_mhz * 1e6)
    return rm * np.power(lamb, 2) + offset


def fit_psi_sky(freqs, S_ratio, pa_fit):
    # psi'(nu) = psi(nu) + psi_sky
    print("Fitting psi_sky")
    Q_askap = np.copy(S_ratio[0])
    popt, pcov = curve_fit(QoverI_askap, pa_fit, Q_askap, p0=[0.95, -0.8])
    L_amp, psi_sky = popt
    L_amp_err, psi_sky_err = np.diag(pcov)

    print(f"L_amp\t= {L_amp:.3f}\t+- {L_amp_err:.5f}")
    print(f"psi_sky\t= {psi_sky:.3f}\t+- {psi_sky_err:.5f}")

    Q_fit = QoverI_askap(pa_fit, L_amp, psi_sky)

    return Q_askap, (L_amp, L_amp_err), (psi_sky, psi_sky_err), Q_fit


def QoverI_askap(pa, L_amp, psi_sky):
    return L_amp * np.cos(2 * pa + psi_sky)


def fit_phi(freqs, S_ratio, pa_fit, L_amp, psi_sky):
    # Solve using matrices
    tan_phi = np.full((freqs.shape[0], 2), np.nan)
    u_prime = S_ratio[1]  # U/I
    v_prime = S_ratio[2]  # V/I
    for i, f in enumerate(freqs.value):
        pa_prime = psi_sky + pa_fit[i]

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

    # fit phi in several overlapping sub-bands, take best fit
    n_bands = 3
    overlap = 2
    band_width = phi.shape[0] // n_bands
    slices = [
        slice(int(i*band_width/overlap), int((i+overlap)*band_width/overlap)) 
        for i in range(n_bands*overlap-(overlap-1))
    ]

    delays = []
    offsets = []
    r2s = []

    for i, s in enumerate(slices):
        pfit, pcov = np.polyfit(freqs[s].value, phi[s], 1, cov=True)
        perr = np.sqrt(np.diag(pcov))

        delay = pfit[0]
        delay_err = perr[0]
        delays += [(delay, delay_err)]

        offset = pfit[1]
        offset_err = perr[1]
        offsets += [(offset, offset_err)]

        p = np.poly1d([delay * (2 * np.pi * 1e6), offset])

        #https://stackoverflow.com/questions/29003241/how-to-quantitatively-measure-goodness-of-fit-in-scipy
        ss_res = np.sum((phi[s] - p(freqs[s].value)**2))
        ss_tot = np.sum((phi[s] - np.mean(phi[s]))**2)
        r2s += [1 - (ss_res / ss_tot)]

        print(f"Band {i+1}: {freqs[s][0]} - {freqs[s][-1]}")
        print(
            f"Delay:  {delay} +- {delay_err} rad us"
        )
        print(
            f"Offset: {offset} +- {offset_err} "
        )

    best_band = np.argmin(r2s)
    print(f"Band {best_band+1} is the best band")

    phi_fit = np.poly1d(
        [delays[best_band][0], offsets[best_band][0]]
    )(freqs.value)

    return phi, delays[best_band], offsets[best_band], phi_fit


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
    print(("rm", rm, "offset", offset))
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
    title=None,
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
        image
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
    ax.set_title(title)
    return ax


def im_plot(
    ax,
    z,
    cbar=False,
    xlabel=None,
    ylabel=None,
    title=None,
    aspect="auto",
    interpolation="none",
    orientation="vertical",
    cmap='viridis',
    **kwargs,
):
    im = ax.imshow(
        z.T[::-1], 
        aspect=aspect, 
        interpolation=interpolation, 
        cmap=cmap,
        **kwargs)
    if cbar:
        plt.colorbar(im, ax=ax, orientation=orientation, label=cbar)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return ax


if __name__ == "__main__":
    _main()
