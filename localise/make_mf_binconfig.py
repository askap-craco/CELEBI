##===============================================##
##===============================================##
## Author: Tyson Dial
## Email: tdial@swin.edu.au
## Created: 11/09/2024 
## Last Updated: 01/10/2024
##
##
## This script takes in High Time resolution data
## from CELEBI and creates the following:
##
## 1. Matched filter binconfig file with onpulse
##    region of burst and off-pulse rfi region
##
## 2. wmask structure that holds {time-dependent,
##    freq-dependent weights, finder and rfi masks}
##
##===============================================##
##===============================================##

# imports
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from math import ceil


class empty:
    pass



##  function to average data    ##
def average(x: np.ndarray, axis: int = 0, N: int = 10, nan = False):

    """
    average in either frequency or time

    Parameters
    ----------
    x: ndarray
       data to average over
    axis: int 
        axis to average over
    N: int
        Averaging/donwsampling factor
    nan : bool, optional
        If True, using nanmean to ignore NaN values in array 'x', by default False

    Returns
    -------
    x: ndarray 
       Averaged data
    
    """

    # if nan if true, will use numpy function that ignores nans in array x
    if nan:
        func = np.nanmean
    else:
        func = np.mean


    if N == 1:
        return x
    

    # either dynamic spectra or time series
    ndims = x.ndim
    if ndims == 1:
        N_new = int(x.size / N) * N
        return func(x[:N_new].reshape(int(N_new / N), N),axis = 1).flatten()
    
    elif ndims == 2:
        if axis == 0:
            #frequency scrunching
            N_new = int(x.shape[0] / N) * N
            return func(x[:N_new].T.reshape(x.shape[1],int(N_new / N), N),axis = 2).T
        
        elif axis == 1 or axis == -1:
            #time scrunching
            N_new = int(x.shape[1] / N) * N
            return func(x[:,:N_new].reshape(x.shape[0],int(N_new / N), N),axis = 2)
        
        else:
            print("axis must be 1[-1] or 0")
            return x
        
    else:
        print("ndims must equal 1 or 2..")
        return x




def get_args():
    """
    Get args
    """

    # arguments
    parser = argparse.ArgumentParser()

    # data arguments
    parser.add_argument('-i', help = "Stokes I ds filepath", type = str) 
    parser.add_argument('-s', help = "Summary filepath with parameters", type = str)
    parser.add_argument('-b', help = "Initial Binconfig file", type = str)
    parser.add_argument('-p', help = "Polyco file", type = str)

    # data processing arguments
    parser.add_argument("--thres", help = "threshold used to define bounds of on-pulse burst based on percentage of maxima", type = float,
                            default = 0.05)
    parser.add_argument("--tN", help = "Factor for averaging in time", type = int, default = 10)
    parser.add_argument("--rms_g", help = "guard between peak and off-pulse region in phase units", type = float, 
                            default = 0.1)
    parser.add_argument("--rms_w", help = "width of off-pulse region in phase units", type = float, default = 0.025)
    parser.add_argument("--rfi_w", help = "Width of rfi region in ms", type = float, default = 4.0)
    parser.add_argument("--rfi_g", help = "Width of guard region between on-pulse and off-pulse in ms", type = float, default = 1.0)

    # weights
    parser.add_argument("--tw", help = "Calculate time-dependent weights to apply", action = "store_true")
    parser.add_argument("--fw", help = "Calculate freq-dependent weights to apply", action = "store_true")

    return parser.parse_args()







def load_files(args):
    """
    Load parameter files and add to args

    Parameters
    ----------
    args : argparse.argumentparser
        arguments for script

    Returns
    -------
    args : argparse.argumentparser
        arguments for script with additional paramters from loaded files
    """

    # load Stokes I data file
    args.ds = np.load(args.i, mmap_mode = 'r')

    # load summary file
    # parameters we want from this file -> 
    parameters = ["FRB name", "cfreq", "bw", "DM_ref_freq", "htr_DM", "corr_DM", "corr_ref_freq", "corr_MJD",
                  "Geocentric delay", "crop_MJD"]
    attr       = ["name", "cfreq", "bw", "DM_ref_freq", "htr_DM", "corr_DM", "corr_ref_freq", "corr_MJD",
                  "geodelay", "crop_MJD"]
    with open(args.s, "r") as file:
        lines = file.readlines()
        first_line, last_line = None, None
        for i, line in enumerate(lines):
            if "GENERAL DATA" in line:
                first_line = i+1    # set first line
            
            if "POSITION" in line:
                last_line = i-1     # set last line
                break
        
        lines = lines[first_line:last_line]

    # split up lines and get parameter names
    params = {sstr.split(':')[0] : sstr.split(':')[1].split()[0] for sstr in lines}

    # get requested parameters
    for i, requested_param in enumerate(parameters):
        if requested_param in params.keys():
            if requested_param == "FRB name":
                setattr(args, attr[i], params[requested_param])
            else:
                setattr(args, attr[i], float(params[requested_param]))
        else:
            print(f"Missing [{requested_param}] in Summary file!")
            sys.exit()
    
    # load binconfig file
    # parameters we want from this file ->
    # [numpolycofiles, polyco_file, scrunchoutput]
    with open(args.b, "r") as file:
        lines = file.readlines()
        for line in lines:
            if "NUM POLYCO FILES" in line:
                args.numpolycofiles = line.split(':')[1].strip()
            if "POLYCO FILE 0" in line:
                args.polycofile = "./craftfrb.polyco"
            if "SCRUNCH OUTPUT" in line:
                args.scrunch_output = line.split(':')[1].strip()

    # load polycofile
    with open(args.p, "r") as file:
        lines = file.readlines()
        pulsar_freq = float(lines[1].split()[1])
        args.pulsar_period = 1/pulsar_freq

    return args






def crop_frb(args):
    """
    Crop FRB using sigma values 

    Parameters
    ----------
    args : argparse.argumentparser
        arguments for script with additional paramters from loaded files

    Returns
    -------
    ds : np.ndarray
        cropped dynamic spectra
    """

    MAX_BINS = 120
    TN_ITERATION = 108
    RMS_THRES = 1.5

    # Find on-pulse and off-pulse regions for FRB and rfi respectivley.

    # first we implement a robust check to see if we are not exceeding a certain number of correlator bins.
    while True:
        # average data with current time resolution
        ds = average(args.ds, axis = 1, N = args.tN)
        ds = average(ds, axis = 0, N = 4, nan = True)

        # get on-pulse region bounds, first find max
        t = np.nanmean(ds, axis = 0)
        t_peak = np.argmax(t)

        # get rms as a secondary threshold
        rms_guard_nsamp = int(args.rms_g * t.size)
        rms_nsamp = int(args.rms_w * t.size)
        t_rms = np.std(t[t_peak - rms_guard_nsamp - rms_nsamp:t_peak - rms_guard_nsamp])

        # set all data in t < 1.5 t_rms to zero
        t[t < RMS_THRES * t_rms] = 0.0

        # filter out bins with signal
        t_ind = np.where(np.abs(t/t[t_peak]) > args.thres)
        burst_start_samp = np.min(t_ind)
        burst_end_samp = np.max(t_ind)
        burst_nsamp = burst_end_samp - burst_start_samp

        # calculate number of bins that will be correlated
        rfi_nsamp = int(max(args.rfi_w * 1000 / args.tN, 0))       # AB, 31 Oct 25
        guard_nsamp = int(args.rfi_g * 1000 / args.tN)

        if (burst_nsamp + 2*rfi_nsamp + 2*guard_nsamp) < MAX_BINS:
            print(f"{args.tN} us bins will be used for matched filter correlation")
            break

        print(f"{args.tN} us bins are too fine with {burst_nsamp} on-pulse bins, {guard_nsamp} guard bins and {rfi_nsamp} rfi bins, trying {args.tN + TN_ITERATION}")
        
        # redo with higher time scrunching, add 50 to args.tN for each iteration
        args.tN += TN_ITERATION


    # finally, crop dynamic spectrum, making sure to index at original sample resolution
    ds = ds[:, burst_start_samp - guard_nsamp - rfi_nsamp:burst_end_samp + guard_nsamp + rfi_nsamp]

    print("Crop found at sample positions:")
    print(f"Starting onpulse sample: {burst_start_samp}")
    print(f"Ending onpulse sample: {burst_end_samp}")
    # print(f"DM sweep sample number: {delDM_sweep_samp}")
    print(f"RFI guard sample number: {guard_nsamp} X 2")
    print(f"RFI region sample number: {rfi_nsamp} X 2")

    # make relative
    args.burst_start_samp = burst_start_samp
    args.burst_end_samp = burst_end_samp
    args.guard_nsamp = guard_nsamp
    args.rfi_nsamp = rfi_nsamp
    args.crop_start_samp = burst_start_samp - guard_nsamp - rfi_nsamp
    args.crop_end_samp = burst_end_samp + guard_nsamp + rfi_nsamp

    # bin indexes
    args.finder_start = args.burst_start_samp - args.crop_start_samp
    args.finder_end = args.burst_end_samp - args.crop_start_samp

    # other params
    args.rms_nsamp = rms_nsamp
    args.t_rms = t_rms

    return ds, args












def make_binconfig(ds, args):
    """
    Make binconfig file with time-dependent weights

    Parameters
    ----------
    ds : np.ndarray
        cropped Dynamic spectrum
    args : argparse.argumentparser
        arguments for script with additional paramters from loaded files

    Returns
    -------
    tw : np.ndarray
        time dependent weights from cropped FRB
    fw : np.ndarray
        freq dependent weights from cropped FRB
    """

    # create weight mask object
    # FORMAT OF WEIGHT OUTPUT FILE #TODO - can edit later to add additional weighting profiles
    # -> T_WEIGHTS
    # -> F_WEIGHTS
    # -> T_WEIGHTS MASK FOR WEIGHTING
    # -> FINDER MASK (HTR mask)
    # -> RFI MASK
    # -> FIELD MASK
    wmask = empty()

    # for diagnostics
    diagouts = empty()
    diagouts.ds = ds.copy()

    # create rfi-subtracted crop of data
    if (args.rfi_nsamp > 0):
        mean_rfi = (np.mean(ds[:,:args.rfi_nsamp], axis = 1) 
                    + np.mean(ds[:,args.finder_end + args.guard_nsamp:], axis = 1))/2
        ds -= mean_rfi[:, None]
    diagouts.ds_rfisub = ds.copy()

    # Get time(freq-)dependent weights
    tw = np.zeros(ds.shape[1])
    fw = np.zeros(ds.shape[0])

    if args.tw:
        tw[args.finder_start:args.finder_end] = np.nanmean(ds[:,args.finder_start:args.finder_end], axis = 0)
    else:
        tw[args.finder_start:args.finder_end] = 1.0

    if args.fw:
        fw = np.mean(ds[:,args.finder_start:args.finder_end], axis = 1)
    else:
        fw = np.ones(ds.shape[0])

    # set any flagged channels to zero weights
    fw[np.isnan(fw)] = 0.0

    # set any negative values to 0.0
    fw[fw < 0.0] = 0.0
    tw[tw < 0.0] = 0.0

    # normalise weights to max
    tw /= np.max(tw)
    fw /= np.max(fw)

    # zero any weights less then 1%
    tw[tw < 0.01] = 0.0

    # add rfi bin weights
    tw[:args.rfi_nsamp] = -1
    tw[args.finder_end + args.guard_nsamp:] = -1

    # add zeroth bin (field bin)
    tw = np.concatenate(([0],tw))

    # add time(freq-)weights to container, also add zeroth bin for 
    wmask.tw = tw.copy()
    wmask.fw = fw.copy()
    
    # create finder, rfi and field masks
    # just make finder bin the ds with padded zeroth bin (NOTE: may be removed later)
    wmask.findermask = np.pad(ds, ((0,0),(1,0)))

    wmask.tmask = np.tile(tw, (fw.size, 1))
    wmask.rfimask = wmask.tmask.copy()
    wmask.tmask[wmask.tmask < 0] = 0.0          # remove rfi bins in tmask
    wmask.rfimask[wmask.rfimask != -1] = 0.0    # remove all but rfi bins

    wmask.fieldmask = np.zeros(wmask.tmask.shape)
    wmask.fieldmask[:,0] = 1.0



    # Now construct the binconfig file, steps to take:
    # 1. Correct MJD timestamping for imaging
    # 2. Construct binconfig file with new weights

    # 1. Correction for MJD timestamp, need to account for geometric delay and differences in de-dispersion
    geo_delay_MJD = args.geodelay / 86400
    DM_delay_MJD = (4149.377593 * args.htr_DM * 
                            (1/args.DM_ref_freq**2 - 1/args.corr_ref_freq**2)) / 86400
    
    args.geodelay_ms = args.geodelay * 1000
    args.DM_delay_ms = DM_delay_MJD * 86400 * 1000

    print(f"DM sweep: {4149.377593 * args.htr_DM * (1/args.DM_ref_freq**2 - 1/args.corr_ref_freq**2)}")
                        
    # get MJD of first sample in cropped dynamic spectra, this will be the first phase bin
    bin0_MJD = args.crop_MJD + args.crop_start_samp * args.tN / 8.64e10
    print(f"MJD before correction: {bin0_MJD}")
    print(f"Time before MJD corr: {(bin0_MJD - args.corr_MJD) * 86400} s")
    bin0_MJD += (geo_delay_MJD - DM_delay_MJD)
    args.bin0_MJD = bin0_MJD 
    print(f"MJD after correction: {bin0_MJD}")
    print(f"Time after MJD corr: {(bin0_MJD - args.corr_MJD) * 86400} s")


    # Build new binconfig file
    nbins = ds.shape[1]
    new_binconfig = "mf.binconfig"
    with open(new_binconfig, "w") as binconfig:

        # Write header infomation
        binconfig.write("NUM POLYCO FILES:".ljust(20) + args.numpolycofiles + "\n")
        binconfig.write("POLYCO FILE 0:".ljust(20) + args.polycofile + "\n")
        binconfig.write("NUM PULSAR BINS:".ljust(20) + str(nbins + 1) + "\n")
        binconfig.write("SCRUNCH OUTPUT:".ljust(20) + args.scrunch_output + "\n")

        # write zeroth bin
        bin0_phase = (bin0_MJD - args.corr_MJD) * 86400 / args.pulsar_period
        binconfig.write("BIN PHASE END 0:".ljust(20) + f"{bin0_phase:.7f}" + "\n")
        binconfig.write("BIN WEIGHT 0:".ljust(20) + "1.0\n")
        
        # write the rest of the bins
        del_phase = (1e-6 * args.tN) / args.pulsar_period
        for i in range(nbins):
            binphase = bin0_phase + (i+1) * del_phase
            if binphase > 1.0:
                print(f"Warning: binphase for bin {i} > 1.0, double check the 'fake' pulsar period in the .polyco file!!!")
            binconfig.write(f"BIN PHASE END {i+1}:".ljust(20) + f"{binphase:.7f}" + "\n")
            binconfig.write(f"BIN WEIGHT {i+1}:".ljust(20) + "1.0\n")

    print(f"Written new binconfig file to [{new_binconfig}]")


    # Build new polyco file
    # load old polyco file
    with open(args.p, "r") as file:
        lines = file.readlines()

    new_polyco = "mf.polyco"
    with open(new_polyco, "w") as polyco:

        # update DM to match HTR
        line1 = lines[0].split()
        lines[0] = f"{line1[0]} {line1[1]} {line1[2]} {line1[3]} {args.htr_DM} {line1[5]} {line1[6]}\n"
        polyco.writelines(lines)
    
    print(f"Written new polyco file to [{new_polyco}]")


    return wmask, args





def diagnostics(ds, args, wmask):
    """
    Diagnostics i.e. plotting, saving data files etc.

    """

    # save single plot of Dynspec and time series plot
    fig, AX = plt.subplot_mosaic([['f', 'cas'],['None', 't']], figsize = (10,10), 
                             gridspec_kw={'height_ratios':[8,1], 'width_ratios':[1,10]})

    # remove unused axis
    AX['None'].remove()

    # dynamic spec
    AX['cas'].get_xaxis().set_visible(False)
    AX['cas'].get_yaxis().set_visible(False)
    AX['cas'].set_title("HTR Weights")
    AX['cas'].imshow(ds, aspect = 'auto', extent = [0, ds.shape[1], args.cfreq - args.bw/2,   
                                                args.cfreq + args.bw/2])

    # time weights
    AX['t'].set_xlabel("Time Bins", fontsize = 16)
    AX['t'].get_yaxis().set_visible(False)
    tw = wmask.tw.copy()
    tw[tw < 0] = 0.0
    AX['t'].imshow(tw[1:].reshape(1, ds.shape[1]), aspect = 'auto',
                        extent = [0, ds.shape[1], 0, 1])

    # freq weights
    AX['f'].get_xaxis().set_visible(False)
    AX['f'].set_ylabel("Freq [MHz]", fontsize = 16)
    AX['f'].imshow(wmask.fw.reshape(ds.shape[0], 1), aspect = 'auto',
                        extent = [0, 1, args.cfreq - args.bw/2, args.cfreq + args.bw/2])

    # adjust figure
    fig.tight_layout()
    fig.subplots_adjust(hspace = 0, wspace = 0)

    # save figure
    plt.savefig("mf_weights.png")

    # make plot of masks
    plt.figure(figsize = (10, 10), layout = "constrained")
    plt.imshow(wmask.tmask + wmask.rfimask, aspect = 'auto', extent = [0, ds.shape[1], args.cfreq - args.bw/2,   
                                                args.cfreq + args.bw/2], cmap = 'twilight_shifted') 
    plt.xlabel("Time Bins", fontsize = 16)
    plt.ylabel("Freq [MHz]", fontsize = 16)
    plt.colorbar()

    plt.savefig("finder_rfi_bins.png")


    # Make figure of time series with rms marker and rfi bounds
    fig3, ax3 = plt.subplots(figsize = (10,10), layout = "constrained")

    tMAX = np.max(np.nanmean(ds, axis = 0))

    ax3.plot(np.nanmean(ds, axis = 0) / tMAX)
    ax3.plot([0, ds.shape[1]-1], [args.t_rms / tMAX]*2, 'k--', label = "rms")
    ylim = ax3.get_ylim()
    xlim = ax3.get_xlim()

    # rfi bin markers
    ax3.plot([1, 1], ylim, 'r--')
    ax3.plot([args.rfi_nsamp]*2, ylim, 'r--', label = "RFI bins")
    ax3.plot([args.finder_end + args.guard_nsamp]*2, ylim, 'r--')
    ax3.plot([args.finder_end + args.guard_nsamp + args.rfi_nsamp - 1]*2, ylim, 'r--')    
    
    # finder bin markers
    ax3.plot([args.finder_start]*2, ylim, 'm--', label = "FINDER bins")
    ax3.plot([args.finder_end]*2, ylim, 'm--')

    # S/N markers
    ax3.plot(xlim, [0.8, 0.8], label = "mf_thres = 0.8", alpha = 0.8, linestyle = ':')
    ax3.plot(xlim, [0.5, 0.5], label = "mf_thres = 0.5", alpha = 0.8, linestyle = ':')
    ax3.plot(xlim, [0.3, 0.3], label = "mf_thres = 0.3", alpha = 0.8, linestyle = ':')
    ax3.plot(xlim, [0.2, 0.2], label = "mf_thres = 0.2", alpha = 0.8, linestyle = ':')
    ax3.plot(xlim, [0.15, 0.15], label = "mf_thres = 0.15", alpha = 0.8, linestyle = ':')
    ax3.plot(xlim, [0.1, 0.1], label = "mf_thres = 0.1", alpha = 0.8, linestyle = ':')
    ax3.plot(xlim, [0.05, 0.05], label = "mf_thres = 0.05", alpha = 0.8, linestyle = ':')
    ax3.plot(xlim, [0.03, 0.03], label = "mf_thres = 0.03", alpha = 0.8, linestyle = ':')

    ax3.legend()

    ax3.set_ylim(ylim)
    ax3.set_xlim(xlim)
    ax3.set_xlabel("Time Bins", fontsize = 16)
    ax3.set_ylabel("(arb.)")

    plt.savefig("htr_t.png")
    

    # save data in wmask class
    with open("mf_wmask.npy", 'wb') as file:

        np.save(file, wmask.tw)             # time-dependent weights
        np.save(file, wmask.fw)             # freq-dependent weights
        np.save(file, wmask.tmask)          # time weight mask
        np.save(file, wmask.findermask)     # finder mask
        np.save(file, wmask.rfimask)        # rfi mask
        np.save(file, wmask.fieldmask)      # field mask
    

    # save info file about new files
    with open("mf_info.txt", "w") as file:
        file.write(f"Matched Filter information for [{args.name}]\n")
        file.write(f"-"*40 + "\n")
        file.write(f"Start MJD of mf binconfig file: {args.bin0_MJD}\n")
        file.write(f"Width of binconfig file: {args.tN * ds.shape[1] * 1e-3} ms\n")
        file.write(f"Number of bins: {ds.shape[1]}\n")
        file.write(f"Number of FINDER bins: {args.finder_end - args.finder_start}\n")
        file.write(f"Number of RFI bins: {2 * args.rfi_nsamp}\n")
        file.write(f"Geometric delay: {args.geodelay_ms} ms\n")
        file.write(f"DM delay (due to difference in DM reference frequency used between correlation and beamforming): {args.DM_delay_ms} ms\n")
        file.write(f"Pulsar period: {args.pulsar_period} s\n")










if __name__ == "__main__":
    # main code block, entry point

    # get args
    args = get_args()

    
    # load files
    args = load_files(args)

    # # crop FRB
    ds, args = crop_frb(args)

    # # make binconfig
    wmask, args = make_binconfig(ds, args)

    # plot/save mf data
    diagnostics(ds, args, wmask)

    # DONE
