##################################################
# Author:   Tyson Dial                           #
# Contributions: Based on work done by Danica    #
# Scott and Apurba Bera                          #
# Email:    tdial@swin.edu.au                    #
# Date (created):     23/04/2024                 #
# Date (updated):     26/04/2024                 #
##################################################
# make HTR mosaic plot                           #          
#                                                #
##################################################
from mmap import mmap
import matplotlib as mpl
mpl.use('agg')

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
from os import path


default_col = plt.rcParams['axes.prop_cycle'].by_key()['color']





def get_args():
    parser = ArgumentParser(description="Plots FRB Stokes dynamic spectra")

    # for loading in the right Dynspec files
    parser.add_argument("-s", type=str, required=True, help="File containing names of dynspec files")

    # for correct frequency ranges
    parser.add_argument("-f", type=float, help="Central frequency (MHz)")

    # for labeling purposes
    parser.add_argument("-l", "--label", type=str, help="FRB label")
    parser.add_argument("-d", "--DM", type=str, help="Dedispersion DM")

    # candidate/time bin files for getting burst MJD
    parser.add_argument("-c", type=str, help="Optimal FRB candidate")
    parser.add_argument("-t", type=float, help="MJD start time")

    # chan flagging
    parser.add_argument("--chanlists", help = "path to dir of files for static channel masking", type = str)

    # mosaic options
    parser.add_argument("--t_panels", help = "Time resolutions to process for HTR mosaic",
                        nargs='+', default = [1, 3, 10, 30, 100, 300, 1000], type = int)

    return parser.parse_args()






def t_average(x, N):
    """
    Average x array by some factor N
    """
    if N == 1:
        return x.copy()

    len_new = (x.shape[1]//N)*N

    return np.mean(x[:,:len_new].reshape(x.shape[0], len_new//N, N), axis = 2)






def load_data(args):
    """
    Get elapsed time in 1us resolution
    """

    # open cand file
    with open(args.c) as cand_file:
        cand = list(map(float, cand_file.readlines()[1].split()))

    # get diff from start MJD to cand MJD
    # NOTE: could improve in the future using beamforming MJD 

    mjd_cand = cand[-2]

    t_elapsed_ms = (mjd_cand - args.t)*86400000

    # convert to ms, take the mid point of the ms bin, i.e. take away 0.5 ms
    t_elapsed_us = (t_elapsed_ms - 0.5) * 1000

    # load in stokes spectra
    stk = {}

    # get filenames of stokes dynspec
    with open(args.s) as fnames_file:
        fnames = fnames_file.readlines()
        fnames = [fname.strip() for fname in fnames]

    # load dynspec as memory map for memory efficiency
    for i, S in enumerate("IQUV"):
        stk[S] = np.load(fnames[i], mmap_mode = "r")

    # convert us time -> bin
    t_burst_bin = int(t_elapsed_us)

    return stk, t_burst_bin








def flag_chan(dsI, on_pulse, flag_thresh, tN):
    """
    zap channels dynamicall based on RMS
    NOTE: Code developed by Apurba Bera (2 Apr 2023) and cleaned by Tyson Dial (26 Apr 2024)
    
    Paramters
    ---------
    stk          -> dictionary of stokes dynspec
    on_pulse     -> Slice object of on-pulse region
    flag_thresh  -> Outlier threshold in units of SD
    tN           -> Averaging factor in time (integer)
    """
    
    # crop out on-pulse region
    dsI[:, on_pulse] = np.nan

    # use stokes I dynamic spectrum
    dsI = t_average(dsI, tN)
    
    # channel mask, this is used to mask out bad channels (False entries)
    chanmask = np.ones(dsI.shape[0], dtype = bool)

    # flag based on noise in each channel
    fI_std = np.nanstd(dsI, axis = 1)
    med_rms = np.nanmedian(fI_std)
    mad_rms = 1.48 * np.nanmedian(np.abs(fI_std - med_rms))
    chan2flag = np.where(fI_std > (med_rms + flag_thresh*mad_rms))[0]
    chanmask[chan2flag] = False

    # instrumentation zapping 
    # get flag file for given bandwidth, NOTE: This could be put in the nextflow script instead??
    askap_badchan_file = path.join(args.chanlists, "htrchanlist_low.txt")
    if args.f > 1100.0:
        askap_badchan_file = path.join(args.chanlists, "htrchanlist_mid.txt")
    if args.f > 1500.0:
        askap_badchan_file = path.join(args.chanlists, "htrchanlist_high.txt")

    # flag bad channels within bandwidth
    askap_chan2flag = np.loadtxt(askap_badchan_file)
    if askap_chan2flag.shape[0] > 2:
        for i in range(2, askap_chan2flag.shape[0]):
            chanmask[int(round(askap_chan2flag[i,0])):int(round(askap_chan2flag[i,1]))+1] = False
    
    return chanmask










def plot_htr(args, stk, t_burst_bin):
    """
    create mosaic of frb burst at different t resolutions for 
    IQUV stokes dynamic spectra

    Parameters
    ----------
    args        - > function arguments
    stk         - > dictionary of stokes dynspec
    tburst_bin  - > finder bin index (where the burst is approximately to 1ms resolution)
    """

    # can change later
    t_arr = args.t_panels
    nsamp = 100


    # construct mosaic figure
    num = len(t_arr)
    pmax = max(t_arr)

    # create figure
    axes_handles = [[f"{S}{t}" for t in t_arr + ["f"]] for S in "tIQUV"]
    x_plot_w = 2*7/num
    fig, AX = plt.subplot_mosaic(axes_handles, figsize = (18,12),
            gridspec_kw = {"height_ratios": [1,2,2,2,2], "width_ratios": [x_plot_w]*num+[1]})

    # channel zap 
    rough_on_pulse = slice(t_burst_bin - int(1.2*1000*nsamp), t_burst_bin + int(1.2*1000*nsamp) + 1)
    chanmask = flag_chan(stk['I'].copy(), rough_on_pulse, 10.0, 1000)

    # preprocess stokes dynspec
    # find robust peak in data
    tI = np.mean(stk['I'][chanmask], axis = 0)
    tI = t_average(tI.reshape(1, tI.size), pmax).flatten()
    peak = np.argmax(tI) * pmax

    # crop data
    stk_data = {}
    on_pulse = slice(peak - int(1.2*pmax*nsamp), peak + int(1.2*pmax*nsamp) + 1)
    for S in "IQUV":
        stk_data[S] = stk[S][:,on_pulse].copy()

    # flag stokes data
    for S in "IQUV":
        stk_data[S][~chanmask] = np.nan

    
    # loop through t factors
    freqs = np.linspace(args.f + 167.5, args.f - 167.5, 336)
    flim = [freqs[0], freqs[-1]]
    S_col = ['']
    for j, t in enumerate(t_arr):

        # loop through stokes IQUV
        for k,S in enumerate("IQUV"):
            # diagnostics
            print(f"Processing: t = {t}, S = {S}: {k+4*j+1}/{4*len(t_arr)}")

            # crop and average by factor t
            ds = t_average(stk_data[S], t)

            # max point w.r.t power spectrum (I)
            if S == "I":
                peak = np.argmax(np.nanmean(ds, axis = 0))

            ds = ds[:,peak - nsamp:peak + nsamp + 1]

            # get time bins of burst in [ms]
            times = np.linspace(0.5*t, t*ds.shape[1] - 0.5*t, ds.shape[1])/1e3
            tlim = [times[0], times[-1]]

            # plot dynamic spectrum
            AX[f"{S}{t}"].imshow(ds, aspect = 'auto', 
                                 extent = [*tlim, *flim[::-1]])
            
            # plot stokes time series
            AX[f"t{t}"].plot(times, np.nanmean(ds, axis = 0), color = default_col[k],
                                label = S)
            

            # plot stokes I freq spectra
            if t == t_arr[-1]:
                # nsamp is our maximum, from the way we cropped the data
                AX[f"{S}f"].plot(ds[:,nsamp], freqs, color = 'k')
                AX[f"{S}f"].plot([0.0]*2, flim, '--k', linewidth = 1.5)
                AX[f"{S}f"].set_ylim(flim[::-1])

            
            # set axes dynsec axis properties
            if t == t_arr[0]:
                AX[f"{S}{t}"].set_ylabel("Frequency [MHz]")
            else:
                AX[f"{S}{t}"].get_yaxis().set_visible(False)
            
            if S == "V":
                AX[f"{S}{t}"].set_xlabel("t offset [ms]")
            else:
                AX[f"{S}{t}"].get_xaxis().set_visible(False)
            

        # set stokes time series axis properties
        AX[f"t{t}"].get_xaxis().set_visible(False)
        AX[f"t{t}"].set_xlim(tlim)
        AX[f"t{t}"].set_title(f"{t} $\\mu$s", fontsize = 20)
        if t == t_arr[0]:
            AX[f"t{t}"].set_ylabel("Flux Density (arb.)")
        else:
            AX[f"t{t}"].get_yaxis().set_visible(False)
        if t == t_arr[-1]:
            AX[f"t{t}"].legend()
            AX[f"t{t}"].get_legend().set(bbox_to_anchor=(1.03, 0.85))


    # set stokes freq spectra axis properties
    for S in "IQUV":
        AX[f"{S}f"].get_xaxis().set_visible(False)
        AX[f"{S}f"].set_yticks([])
        AX[f"{S}f"].yaxis.set_label_position("right")
        AX[f"{S}f"].set_ylabel(S, rotation = 0, fontsize = 20, labelpad = 20)

    # re-calibrate ylim of time series plots on top of mosaic plot
    ylim = AX[f"t{t_arr[-1]}"].get_ylim()
    for t in t_arr:
        AX[f"t{t}"].set_ylim(ylim)

    # remove corner axes
    AX[f"tf"].set_axis_off()

    # final figure adjustments
    fig.tight_layout()
    fig.subplots_adjust(hspace = 0, wspace = 0)  


    # save figure
    plot_filename = f"{args.label}_IQUV_{args.DM}.png"
    print(f"Saving figure as: {plot_filename}")
    plt.savefig(plot_filename)


    # save flagged channels to file
    np.savetxt(f"{args.label}_channel_mask.txt", chanmask, fmt='%d')


    # saved cropped data to file
    for S in "IQUV":
        print(f"Saving crop of Stokes {S} dynspec as: {args.label}_{args.DM}_ds{S}_crop.npy")
        np.save(f"crops/{args.label}_{args.DM}_ds{S}_crop.npy", stk_data[S])
    








if __name__ == "__main__":
    # main block of code

    # get arguments
    args = get_args()

    # get burst starting time and load in stokes dynspecs
    stk, t_burst_bin = load_data(args)

    # plot data
    plot_htr(args, stk, t_burst_bin)
