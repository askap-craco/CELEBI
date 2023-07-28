from mmap import mmap
import matplotlib as mpl
mpl.use('agg')

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np


def _main():
    args = get_args()

    with open(args.s) as fnames_file:
        fnames = fnames_file.readlines()
        fnames = [fname.strip() for fname in fnames]
    
    with open(args.c) as cand_file:
        cand = list(map(float, cand_file.readlines()[1].split(" ")))
    
    plot(args, fnames, cand)


def get_args():
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
    parser.add_argument(
        "-x", type=str, help="X fine spectrum"
    )
    parser.add_argument(
        "-y", type=str, help="Y fine spectrum"
    )
    parser.add_argument(
        "-c", type=str, help="Optimal FRB candidate"
    )
    parser.add_argument(
        "-t", type=str, help="Time axis in MJD (1 ms steps)"
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
    init_t,
    facs=[1, 3, 10, 30, 100, 300],
    time_range=100,
    labels=None,
    fname=None,
):
    fig = plt.figure(figsize=(2.5 * len(facs), 2.5 * len(ds_list)))
    spec = gs.GridSpec(nrows=len(ds_list) + 1, ncols=3 * len(facs) + 1)

    abs_peak = None
    peaks = []
    windows = []

    # find peak within init_range, which is centred on expected FRB MJD
    init_range = slice(init_t-time_range, init_t+time_range, 1)

    prof_axs = [fig.add_subplot(spec[0, 3 * j : 3 * (j + 1)]) for j in range(len(facs))]

    for i, ds in enumerate(ds_list):
        print(f"{i}")
        j = len(prof_axs) - 1
        for k, dt in enumerate(facs[::-1]):
            print(f"  {j}, {dt}")
            if i == 0:  # Stokes I, identify peak
                ds_red = reduce(ds, dt, axis=1)
                ts_red = np.sum(ds_red, axis=0)
                if abs_peak is None:
                    first_peak = np.argmax(ts_red)#[init_range])+init_range.start
                    abs_peak = (first_peak * dt,
                                first_peak * dt+dt+1)
            
                peaks.append(
                    abs_peak[0] // dt + 
                    np.argmax(ts_red[abs_peak[0]//dt:abs_peak[1]//dt])
                )

                windows.append(slice(
                    max(0, peaks[k] - time_range)*dt,
                    min(ds_red.shape[1], peaks[k] + time_range)*dt,
                    1
                ))
                peak = peaks[k]
                window = windows[k]
                plot_ds = ds_red[
                    :,
                    max(0, peak - time_range) : min(ds_red.shape[1], peak + time_range)
                ]
                plot_ts = np.sum(plot_ds, axis=0)

            else:   # Q, U, V
                peak = peaks[k]
                window = windows[k]
                plot_ds = reduce(ds[:, window], dt, axis=1)
                plot_ts = np.sum(plot_ds, axis=0)

            print(f"  peak = {peak}")

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

            if k == 0:  # lowest time resolution
                if i == 0:  # Stokes I
                    plot_spec = ds_red[:, peak]
                else:
                    plot_spec = plot_ds[:, time_range]

            j -= 1
            prev_fac = dt

        # spectrum column
        freqs = np.linspace(f0 - 336 / 2, f0 + 336 / 2, plot_spec.shape[0])

        spec_ax = fig.add_subplot(spec[i + 1, -1])
        spec_ax.step(plot_spec, freqs[::-1], c="k", lw=1)
        spec_ax.axvline(0, c="k", lw=1, ls="--")

        spec_ax.set_yticks([])
        spec_ax.set_xticks([])

        spec_ax.set_ylim(freqs[0], freqs[-1])

        spec_ax.yaxis.set_label_position("right")
        spec_ax.set_ylabel(labels[i], rotation=0, fontsize=20, labelpad=20)

    prof_axs[-1].legend(bbox_to_anchor=(1, 0.5), loc="center left")

    plt.subplots_adjust(
        wspace=0, hspace=0, left=0.075, right=0.95, top=0.975, bottom=0.05
    )

    if fname:
        plt.savefig(fname, dpi=200.0)

    print("Done with dynamic spectra")

    return peaks[-1]  # 1 us peak


def crop(stks, peak, dt, labels, time_range, prefix, c=None):
    if c is None:
        if time_range is None:
            c = slice(0, -1, 1)
        else:
            c = slice(
                max(0, (peak - time_range//2)//dt),
                min((peak + time_range//2)//dt, stks[0].shape[1]//dt), 
                1
            )

    for j, s in enumerate(stks):
        s_red = reduce(s, dt, axis=1) if dt > 1 else s
        # .T twice so this works on 1D arrays too
        np.save(f"{prefix}{labels[j]}.npy", s_red.T[c].T)
    
    return c


def plot_ts(ds, idx_1ms, fname):
    ts = np.sum(ds, axis=0)
    ts = reduce(ts, 1000)
    plt.figure(figsize=(10, 5))
    plt.plot(ts)
    plt.axvline(idx_1ms, c='r')
    plt.xlabel("Time (ms)")
    plt.ylabel("I")
    plt.tight_layout()
    plt.savefig(fname)


def plot_stokes_mask_chan(dsall, idx_1ms, outhresh, extims, chanlist, fname, fpltname):
    """
    Function to mask noisy channels and plot 'clean' Stokes parameters
    
    Inputs are -- All Stokes dynamic spectrum
                  
                  location of the burst
                  
                  Outlier threshold in units of SD
                  
                  Time to exclude around burst in units of ms
                  
                  File containing known bad channels
                  
                  File name for masked channels
                  
                  File name for output plot
                                                       AB 27 Feb 2023
                                              Revised  AB  2 Apr 2023
    """

    dspec1ms  = reduce(dsall[0], 1000, axis=1)
    qdspec1ms = reduce(dsall[1], 1000, axis=1)
    udspec1ms = reduce(dsall[2], 1000, axis=1)
    vdspec1ms = reduce(dsall[3], 1000, axis=1)
    dspecdumm = np.copy(dspec1ms)
    dspecdumm[:,int(round(idx_1ms-extims)):int(round(idx_1ms+extims))] = np.nan

    nchans      =  dspec1ms.shape[0]
    chanmask    =  np.ones(nchans, dtype=int)
    chanstoflag =  np.loadtxt(chanlist)
    if(chanstoflag.shape[0]>2):
        for i in range(2,chanstoflag.shape[0]):
            chanmask[int(round(chanstoflag[i,0])):int(round(chanstoflag[i,1]))+1]=0

    for cc in range(0,nchans):
        if(chanmask[cc]==0):
            dspecdumm[cc] = np.nan
            dspec1ms[cc]  = np.nan
            qdspec1ms[cc]  = np.nan
            udspec1ms[cc]  = np.nan
            vdspec1ms[cc]  = np.nan

    chanrmsarr  = np.nanstd(dspecdumm, axis=1)
    medrms      = np.nanmedian(chanrmsarr)
    madrms      = 1.48*np.nanmedian(np.abs(chanrmsarr - medrms))
    badchans    = np.where(chanrmsarr > (medrms + outhresh*madrms))[0]
    chanmask[badchans] = 0
    np.savetxt(fname, chanmask, fmt='%d')

    dspec1ms[badchans] = np.nan
    dspecdumm[badchans]= np.nan
    ts1ms = np.nanmean(dspec1ms,axis=0)
    tsdumm= np.nanmean(dspecdumm,axis=0)
    tsrms = np.nanstd(tsdumm)
    tspeak= np.nanmax(ts1ms)
    peaksn= tspeak/tsrms

    ldspec1ms = np.sqrt(qdspec1ms**2 + udspec1ms**2)
    ldspec1ms[badchans] = np.nan
    vdspec1ms[badchans] = np.nan

    qts1ms = np.nanmean(qdspec1ms,axis=0)
    uts1ms = np.nanmean(udspec1ms,axis=0)
    lts1ms = np.nanmean(ldspec1ms,axis=0)
    ltsdummy = np.copy(lts1ms)
    ltsdummy[int(round(idx_1ms-extims)):int(round(idx_1ms+extims))] = np.nan
    lts1ms = lts1ms - np.nanmean(ltsdummy)
    vts1ms = np.nanmean(vdspec1ms,axis=0)
    
    pats1ms = 0.5*np.arctan(uts1ms/qts1ms)

    fig=plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(2,2,3)
    plt.imshow(dspec1ms,origin='lower',aspect='auto',interpolation='none',vmin=-2.0*medrms,vmax=3.0*medrms)
    plt.xlim([int(round(idx_1ms-52)),int(round(idx_1ms+52))])
    plt.ylabel('Channel')

    ax2 = fig.add_subplot(2,2,1)
    plt.plot(ts1ms, 'b-')
    plt.text(0.05*len(ts1ms), 0.8*np.amax(ts1ms), "Peak S/N = %.2f"%peaksn)
    plt.xlim([int(round(idx_1ms-52)),int(round(idx_1ms+52))])
    plt.xlabel("Time (ms)")
    plt.ylabel("Stokes I")

    ax3 = fig.add_subplot(2,2,2)
    plt.plot(ts1ms, 'k-', label='I')
    plt.plot(lts1ms, 'b-', label='L')
    plt.plot(vts1ms, 'r-', label='v')
    #plt.text(0.05*len(ts1ms), 0.8*np.amax(ts1ms), "Peak S/N = %.2f"%peaksn)
    plt.xlim([int(round(idx_1ms-22)),int(round(idx_1ms+22))])
    plt.xlabel("Time (ms)")
    plt.ylabel("I, L, V")
    plt.legend(loc="upper right")

    ax4 = fig.add_subplot(2,2,4)
    plt.plot(pats1ms*180.0/np.pi, 'b*')
    #plt.text(0.05*len(ts1ms), 0.8*np.amax(ts1ms), "Peak S/N = %.2f"%peaksn)
    plt.xlim([int(round(idx_1ms-22)),int(round(idx_1ms+22))])
    #plt.xlabel("Time (ms)")
    plt.ylabel("PA (deg)")

    plt.tight_layout()
    plt.savefig(fpltname)
    plt.close()

    return 0



def plot(args, stokes_fnames, cand):
    stks = [np.load(f, mmap_mode="r") for f in stokes_fnames]

    # start crops at expected MJD of FRB - identify index to use
    facs = [1, 3, 10, 30, 100, 300, 1000] # reduction factors --> dt in us

    t_mjd = np.load(args.t)
    mjd_cand = cand[-2]
    idx_1ms = np.argmin(np.abs(t_mjd - mjd_cand))

    # IQUV over four timescales
    peak = plot_IQUV_dts(
        stks, 
        args.f,
        idx_1ms,
        facs=facs,
        labels=["I", "Q", "U", "V"], 
        fname=f"{args.label}_IQUV_{args.DM}.png"
    )

    # Full series at 1 ms
    c_full = crop(
        stks,
        peak,
        1000,
        ["I", "Q", "U", "V"],
        None,    # Full range
        f"crops/{args.label}_{args.DM}_1ms_"
    )

    # +- 50 ms at 50 us
    c_100ms_50us = crop(
        stks,
        peak,
        50,
        ["I", "Q", "U", "V"],
        100*1000,
        f"crops/{args.label}_{args.DM}_50us_"
    )
    with open("50us_crop_start_s.txt", "w") as f:
        f.write(f"{c_100ms_50us.start*50e-6}")

    # +- 10 ms at 1 us
    c_20ms_1us = crop(
        stks,
        peak,
        1,
        ["I", "Q", "U", "V"],
        20*1000,
        f"crops/{args.label}_{args.DM}_1us_"
    )

    # Crop 1D XY time series into +- 5 ms at 3 ns
    X = np.load(args.x, mmap_mode="r")
    Y = np.load(args.y, mmap_mode="r")
    c_1D = slice(
        (c_20ms_1us.start+5000)*336,
        (c_20ms_1us.stop-5000)*336,
        1
    )
    crop(
        [X, Y],
        peak*336,
        1,
        ["X", "Y"],
        None,
        f"crops/{args.label}_{args.DM}_",
        c=c_1D
    )

    # plot full Stokes I time series at 1 ms time resolution
    plot_ts(stks[0], idx_1ms, f"{args.label}_I_{args.DM}.png")

    chanlist = '/fred/oz002/askap/craft/craco/CELEBI/flagging/htrchanlist_low.txt'
    if(args.f > 1100.0):
        chanlist = '/fred/oz002/askap/craft/craco/CELEBI/flagging/htrchanlist_mid.txt'
    if(args.f > 1500.0):
        chanlist = '/fred/oz002/askap/craft/craco/CELEBI/flagging/htrchanlist_high.txt'
    # Mask bad channels and calculate S/N
    plot_stokes_mask_chan(stks, idx_1ms, 10.0, 50.0, \
                    chanlist, \
                    f"{args.label}_channel_mask.txt", \
                    f"{args.label}_I_{args.DM}_masked.png")


if __name__ == "__main__":
    _main()
