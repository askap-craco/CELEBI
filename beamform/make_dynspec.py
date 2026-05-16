##################################################
# Author:   Tyson Dial                           #
# Email:    tdial@swin.edu.au                    #
# Date (created):     15/01/2024                 #
# Date (updated):     16/01/2024                 #
##################################################
# Make Dynspecs of stokes IQUV [make_dynspecs.py]#          
#                                                #
# This script makes dynamic spectra of stokes    #
# IQUV with baseline corrections.                #
##################################################

## Imports 
import numpy as np
from copy import deepcopy
from scipy.fft import fft
from math import ceil
import matplotlib.pyplot as plt

## import basic libraries
import argparse, sys
from os import path, mkdir
import shutil


# stokes functions
def stk_I(x,y):
    """
    Make stokes I
    """
    if x is None: return np.abs(y)**2
    if y is None: return np.abs(x)**2
    return np.abs(x)**2 + np.abs(y)**2

def stk_Q(x,y):
    """
    Make stokes Q
    Stokes Q is negated due to nature of ASKAP PAF, 
    this just means switching x and y around.
    """
    
    return np.abs(y)**2 - np.abs(x)**2

def stk_U(x,y):
    """
    Make stokes U
    """

    return 2 * np.real(np.conj(x)*y)

def stk_V(x,y):
    """
    Make stokes V
    """

    return 2 * np.imag(np.conj(x)*y)







## additional functions
def average(x, axis = 0, N = 1, nan = False):
    """
    scrunch data in either f or t

    ##==== inputs ====##
    x:              data
    N:              scrunch factor
    axis:           axis to scrunch in
    nan:            if nan in inputs

    ##==== outputs ====## 
    x:              scrunched data

    """

    if nan:
        avg_func = np.nanmean
    else:
        avg_func = np.mean

    if axis == -1 or axis == 1:
        # time
        t_new = (x.shape[-1]//N) * N
        shape_new = list(x.shape)
        shape_new[-1] = t_new // N
        shape_new.append(N)

        return np.mean(x[...,:t_new].reshape(shape_new), axis = -1)
    
    elif axis == 0:
        # frequency
        f_new = (x.shape[0]//N) * N
        shape_new = list(x.shape)[::-1]
        shape_new[-1] = f_new // N
        shape_new.append(N)

        return np.mean(x.T[...,:f_new].reshape(shape_new), axis = -1).T

    else:
        print("invalid scrunching axis")
        return x
















def get_args():
    """
    Info:
        Get arguments passed during script call


    Args:
        args: Arguments for POLCAL.py script

    """

    parser = argparse.ArgumentParser(
        description = "Fit for pol cal solutions"
    )

    ## data arguments
    parser.add_argument("-x", help = "X pol time series", type = str, default = None)
    parser.add_argument("-y", help = "Y pol time series", type = str, default = None)
    parser.add_argument("--pols", nargs='+', required=True, help="List of expected polarisations (e.g., X Y or just X)")
    parser.add_argument("--nFFT", help = "Number of frequency channels for final dynspec", 
                        type = int, default = 336)
    parser.add_argument("--bline", help = "Apply baseline correction", action = "store_true")

    # chan flagging
    parser.add_argument("--chanlists", help = "path to dir of files for static channel masking", type = str)
    parser.add_argument("--do_chanflag", help = "Do channel flagging based on channel noise", action = "store_true")


    ## data reduction arguments
    parser.add_argument("--sigma", help = "S/N threshold for baseline correction", type = float, default = 5.0)
    parser.add_argument("--baseline", help = "Width of rms crops in [ms]", type = float, default = 50.0)
    parser.add_argument("--tN", help = "Time averaging factor, helps with S/N calculation", type = int, default = 50)
    parser.add_argument("--guard", help = "Time between rms crops and burst in [ms]",
                        type = float, default = 1.0)


    ## Pulsar arguments (Polarisation calibration)
    parser.add_argument("--pulsar", help = "Is HTR products of a pulsar", action = "store_true")
    parser.add_argument("--MJD0", help = "Initial Epoch MJD", type = float, default = None)
    parser.add_argument("--MJD1", help = "Observation MJD", type = float, default = None)
    parser.add_argument("--F0", help = "Initial Epoch pulsar frequency", type = float, default = None)
    parser.add_argument("--F1", help = "Spin-down rate", type = float, default = None)
    parser.add_argument("--DM", help = "Dispersion Measure of Pulsar", type = float, default = None)
    parser.add_argument("--cfreq", help = "Central Frequency", type = float, default = 1271.5)
    parser.add_argument("--bw", help = "bandwidth", type = float, default = 336.0)


    ## output arguments
    parser.add_argument("--ofile", help = "Name of new dynamic spectra", type = str)

    args = parser.parse_args()

    return args

















def load_data(args):
    """
    Info:
        Load in Stokes I, Q, U & V data along with 
        estimates on L/I and V/I coefficients.

    Args
        args: Arguments of POLCAL.py script


    Return
        stk (dict): dict of stokes dynspecs [IQUV]
    
    """

    ## load in stokes data
    expected_pols = args.pols
    pol = {'X': None, 'Y': None}

    if 'X' in expected_pols:
        if not args.x or not path.exists(args.x):
            raise FileNotFoundError("CRITICAL ERROR: 'X' polarisation expected based on --pols, but file is missing or not provided.")
        pol['X'] = np.load(args.x, mmap_mode = 'r')
        
    if 'Y' in expected_pols:
        if not args.y or not path.exists(args.y):
            raise FileNotFoundError("CRITICAL ERROR: 'Y' polarisation expected based on --pols, but file is missing or not provided.")
        pol['Y'] = np.load(args.y, mmap_mode = 'r')


    return pol







def make_ds(xpol, ypol, S = "I", nFFT = 336):
    """
    Info:
        Make Stokes Dynamic spectra with specified stft length 

    Args:
        xpol (mmap): X polarisation
        ypol (mmap): Y polarisation
        stokes (str): Stokes dynspec to make
        nFFT (int): Number of channels 

    Returns:
        ds (ndarray): Raw Dynamic spectrum

    """
    prog_str = f"[Stokes] = {S} with [nFFT] = {nFFT}"

    Stk_Func = {"I":stk_I, "Q":stk_Q, "U":stk_U, "V":stk_V}

    # pre-processing for iterative 
    BLOCK_SIZE = 200e6 # block size in B
    BIT_SIZE = 8       # Bit size in B

    # First need to chop data so that an integer number of FFT windows can
    # be performed. Afterwards, this data will be split into coarse BLOCKS
    # with specific memory constraints. 

    # define parameters
    ref_pol = xpol if xpol is not None else ypol
    nsamps  = ref_pol.size               # original number of samples in loaded dataset
    fnsamps = (nsamps // nFFT) * nFFT    # number of samples after chopping 
    nwind   = fnsamps // nFFT            # number of fft windows along time series

    # memeory block paramters
    nwinb = int(BLOCK_SIZE // (nFFT * BIT_SIZE))    # num windows in BLOCK
    nsinb = int(nwinb * nFFT)                       # num samples in BLOCK
    nblock= int(nwind // nwinb)                     # num BLOCKS

    # create empty array
    ds = np.zeros((nFFT, nwind), dtype = np.float32)

    b_arr = np.empty((0,2), dtype = int)
    i = -1
    for i in range(nblock):
        b_arr = np.append(b_arr, [[i*nwinb,(i+1)*nwinb]], axis = 0)
    # append extra block at end
    if nblock*nsinb < nsamps:
        b_arr = np.append(b_arr, [[(i+1)*nwinb,nwind]], axis = 0)

    # loop over blocks
    for i, b in enumerate(b_arr): # b is bounds of block in nFFT windows
        sb = b * nFFT
        wind_w = b[1] - b[0]

        xp = fft(xpol[sb[0]:sb[1]].copy().reshape(wind_w, nFFT), axis=1) if xpol is not None else None
        yp = fft(ypol[sb[0]:sb[1]].copy().reshape(wind_w, nFFT), axis=1) if ypol is not None else None
        
        ds[:,b[0]:b[1]] = Stk_Func[S](xp, yp).T
        
        # print progress
        print(f"[MAKING DYNSPEC]:    [Progress] = {(i+1)/(nblock+1)*100:3.3f}%:    " + prog_str,
              end = '         \r')

    print("[MAKING DYNSPEC]:    [Progress] = 100.00%:    " + prog_str + "        \n")
    print(f"Made Dynamic spectra with shape [{ds.shape[0]}, {ds.shape[1]}]")

    return ds







def pulse_fold(ds, DM, cfreq, bw, MJD0, MJD1, F0, F1, chanflag, sphase = None, ):
    """
    Info:
        Takes Pulsar dynamic spectrum and folds it, removes periods
        at far left and right sides to avoid band artifacts produced during
        de-dispersion.

    Args:
        ds (ndarray): dynamic spectrum
        MJD0 (float): Initial Epoch MJD
        MJD1 (float): Observation MJD
        F0 (float): initial Epoch Frequency period
        F1 (float): Spin-down rate
        sphase (float): Starting phase of folding, if not given
                        will be estimated (best done using "I" ds)
        chanflag: ndarray
            channels to flag

    Returns:
        ds_f (ndarray): Raw folded Dynamic spectra
        sphase (float): Starting phase of folding from original dynspec

    """
    print("Pulse Folding Dynspec...")

    ## Calculate Period T in [s]
    T = 1/(F0 - F1 * (MJD1 - MJD0)*86400)
    print(f"with period T = {T}")
    dt = 1e-6 * (ds.shape[0]/336) # get time resolution of dynspec

    ## Fold dynamic spectra
    fold_w = int(T / dt)          # fold width in samples (assumed dt = 1 us)
    fold_n_init = int(ds.shape[1]/fold_w)     # initial number of folds

    # get dispersion sweep, calculate number of "broken" pulse periods
    # due to dispersion.
    top_band = args.cfreq + bw/2
    bot_band = args.cfreq - bw/2
    DM_sweep = 4.14938e3 * DM * (1/bot_band**2 - 1/top_band**2) # DM sweep in seconds
    P_sweep = int(DM_sweep/T) + 1
    print(f"DM sweep: {DM_sweep} [ms]")
    print(f"Culling {P_sweep} Periods to the left due to DM sweep")

    if sphase is None:
        # do first run of folding and find maximum, take as starting phase
        fold_n = int((ds.shape[1] - P_sweep * fold_w)/fold_w)
        ts = np.mean(ds[~chanflag, P_sweep * fold_w : (P_sweep + fold_n) * fold_w], axis = 0)
        ts_f = np.mean(ts.reshape(fold_n, fold_w), axis = 0)

        # find phase of maximum given this is enough signal
        phase_max = np.argmax(ts_f)/ts_f.size
        phase_diff = (phase_max - 0.5)
        sphase = int((P_sweep + phase_diff) * fold_w)

    else:
        # put in sample units
        sphase = int(sphase*ds.shape[1])


    # calculate number of folds
    fold_n = int((ds.shape[1]-(sphase+1))/fold_w)     # number of folds
    print(f"Folding {fold_n}/{fold_n_init} periods")


    # reshape to average folds together
    # ignore side periods due to dedispersing
    ds_r = ds[:,sphase:sphase + fold_w * (fold_n)].copy()
    ds_f = np.mean(ds_r.reshape(ds_r.shape[0], (fold_n), fold_w), axis = 1)


    
    return ds_f, sphase / ds.shape[1]







def flag_chan(ds, flag_thresh, tN, args, rbounds = None):
    """
    Flag channels, this algorithm is not perfect since I'm estimating the RFI across the whole buffer,
    but it should give resonable enough results to find the peak and bounds of the burst.

    NOTE: This code is created by [Apurba Bera] and cleaned up by [Tyson Dial]

    Parameters
    ----------
    ds : ndarray
        Dynamic spectrum stokes I
    flag_thresh : float
        channel flagging threshold
    args : _parse_args_
        Aguments passed to script

    Returns
    -------
    flag_ind : ndarray
        channel indicies flagged
    """

    # create boolean array
    chanmask = np.ones(ds.shape[0], dtype = bool)
    
    
    # Channel flagging of known bad ASKAP channels
    if args.chanlists is not None:
        print("Flagging known bad channels on ASKAP...")

        askap_badchan_file = path.join(args.chanlists, "htrchanlist_low.txt")
        if args.cfreq > 1100.0:
            askap_badchan_file = path.join(args.chanlists, "htrchanlist_mid.txt")
        if args.cfreq > 1500.0:
            askap_badchan_file = path.join(args.chanlists, "htrchanlist_high.txt")

        # flag bad channels within bandwidth
        askap_chan2flag = np.loadtxt(askap_badchan_file)
        if askap_chan2flag.shape[0] > 2:
            for i in range(2, askap_chan2flag.shape[0]):
                chanmask[int(round(askap_chan2flag[i,0])):int(round(askap_chan2flag[i,1]))+1] = False

    chanmask_known = chanmask.copy()
    
    # now do auto channel flagging
    if not args.do_chanflag:
        return ~chanmask, ~chanmask_known


    # remove on pulse
    if rbounds is not None:
        # it is assumed that if rbounds is given, proper baseline correction has already been done
        # flagging the on pulse region for better flagging
        onpulse = ds[:, rbounds[0]:rbounds[1]].copy()
        ds[:, rbounds[0]:rbounds[1]] = np.nan
        
        # average
        ds_avg = average(ds, axis = 1, N = tN, nan = True)

        # put onpulse data back, this is done to reduce memory usage
        ds[:, rbounds[0]:rbounds[1]] = onpulse

        # removing the on pulse region
        ds_avg = ds_avg[:, ~np.isnan(ds_avg[0])]

    else:
        
        # rough baseline correction to help channel flagging
        rms_mean = np.nanmedian(ds, axis = 1)[:, None]
        rms_std = np.nanstd(ds, axis = 1)[:, None]
        ds_avg = (ds - rms_mean) / rms_std

        # average
        ds_avg = average(ds_avg, axis = 1, N = tN)


    # calculate channel rms across buffer
    f_std = np.nanstd(ds_avg, axis = 1)
    med_rms = np.nanmedian(f_std)
    mad_rms = 1.48 * np.nanmedian(np.abs(f_std - med_rms))

    # create boolean array
    chanmask = np.ones(ds_avg.shape[0], dtype = bool)
    
    chan2flag = np.where(f_std > (med_rms + flag_thresh*mad_rms))[0]
    chanmask[chan2flag] = False

    

    return ~chanmask, ~chanmask_known





def plot_failed_bline(ds, tN, chanflag):
    """
    Plot dynspec after baseline corrections failed

    Parameters
    ----------
    ds : ndarray
        Dynamic spectra
    """

    fig, ax = plt.subplots(2, 1, figsize = (10,10), gridspec_kw = {'height_ratios':[1,3]})
    ax = ax.flatten()

    ds[chanflag] = np.nan

    ax[0].plot(np.linspace(0, tN*(1e-3)*ds.shape[1], ds.shape[1]), 
                np.nanmean(ds, axis = 0), 'k')
    ax[0].set_ylabel("Flux Density (arb.)")
    ax[0].get_xaxis().set_visible(False)
    
    ax[1].imshow(ds, aspect = 'auto', extent = [0, tN*(1e-3)*ds.shape[1], 0, 336])
    ax[1].set_ylabel("norm bandwidth [MHz]")
    ax[1].set_xlabel("Time [ms]")

    ax[0].set_xlim([0, tN*(1e-3)*ds.shape[1]])

    # final figure adjustments
    fig.tight_layout()
    fig.subplots_adjust(hspace = 0)  

    plt.savefig("fail_bline.png")
    np.save("fail_bline.npy", ds)

    






def baseline_correction(ds, sigma: float = 5.0, guard: float = 1.0, 
                        baseline: float = 50.0, tN: int = 50, chanflag = None, rbounds = None):
    """
    Info:
        Get baseline corrections to the Dynamic spectra data

    Args:
        ds (ndarray): Dynamic spectra
        sigma (float): S/N threshold for bounds
        guard (float): Time in [ms] between rough bounds and rms crop region
        baseline (float): Width of buffer in [ms] to estimate baseline
                          correction
        tN (int): Time Averaging factor for Dynamic spectrum, helps with
                    S/N calculation.
        rbounds (list): Bounds of FRB burst, if Unspecified, the code will do a rough S/N
                        calculation to determine a bursts bounds

    Returns: 
        bs_mean (ndarray): Baseline mean
        bs_std (ndarray): Baseline std
        rbounds (ndarray): Bounds of FRB burst in Phase units

    """      

    print("Applying baseline correction...")

    # static parameters
    rmsg = 0.5   # rms guard in phase difference from peak of burst

    ## calculate time resolution
    dt = 1e-3 * (ds.shape[0]/336) 

    ## ms -> ds time bin converters
    get_units_avg = lambda t : int(ceil(t/(dt * tN)))
    get_units = lambda t : int(ceil(t/dt))

    ## find burst
    if rbounds is None:
        ## Rough normalize 
        ds_r = average(ds, axis = 1, N = tN)
        rmean = np.nanmedian(ds_r, axis = 1)
        rstd = np.nanstd(ds_r, axis = 1)

        ds_rn = ds_r - rmean[:, None]
        ds_rn /= rstd[:, None]

        
        ## find burst bounds
        print("Looking for bounds of burst...")
        # get peak, crop rms and do rough S/N calculation
        t_rn = np.nanmean(ds_rn[~chanflag], axis = 0)
        peak = np.argmax(t_rn)
        rms_w = get_units_avg(baseline)
        rms_crop = np.roll(t_rn, int(rmsg * ds_rn.shape[1]))[peak-rms_w:peak+rms_w]
        rms = np.nanmean(rms_crop**2)**0.5

        # calculate S/N
        t_sn = t_rn / rms
        
        # check if peak in burst is found
        if np.argwhere(t_sn >= sigma).size == 0:
            print("Couldn't find any signal!!! - Check fail_bline.png and fail_bline.npy")
            plot_failed_bline(ds_rn, tN, chanflag)
            return (None,)*3

        rbounds = np.argwhere(t_sn >= sigma)[[0,-1]]/t_sn.size
        rbounds = np.asarray((rbounds*ds.shape[1]), dtype = int)[:,0]

    ## calculate baseline corrections
    
    guard_w = get_units(guard)
    rms_w = get_units(baseline)
    lhs_crop = ds[:,rbounds[0]-guard_w-rms_w:rbounds[0]-guard_w]
    rhs_crop = ds[:,rbounds[1]+guard_w:rbounds[1]+guard_w+rms_w]
    bl_crop = np.concatenate((lhs_crop, rhs_crop), axis = 1)


    bs_mean = np.mean(bl_crop, axis = 1)
    bs_std = np.std(bl_crop, axis = 1)


    return bs_mean, bs_std, rbounds






def plot_bline_diagnostic(ds, rbounds, chanflag, args, label = ""):
    """
    Generate plot of baseline correction performed

    """


    # create figure and axes
    fig, AX = plt.subplots(2, 1, figsize = (8, 12))
    AX = AX.flatten()
    
    ## calculate time resolution
    dt = 1e-3 * (ds.shape[0]/336) 

    ## ms/ or 1000 x dt -> ds time bin converter
    get_units = lambda t : int(ceil(t/dt))

    # get full rbounds and baseline crop as well as a bit of leg room
    crop_start = rbounds[0] - get_units(args.guard + 1.2*args.baseline)
    crop_end = rbounds[1] + get_units(args.guard + 1.2*args.baseline)


    # crop
    ds_crop = average(ds[:,crop_start:crop_end], axis = 1, N = args.tN)

    # flag channels
    ds_crop[chanflag] = np.nan

    # get time series
    t_crop = np.nanmean(ds_crop, axis = 0)

    # get time axis in ms/ or 1000 x dt
    x_crop = np.linspace(0, dt*args.tN*t_crop.size, t_crop.size)


    ## plot
    AX[0].plot(x_crop, t_crop, color = 'k')
    ylim = AX[0].get_ylim()
    AX[0].plot([0.2*args.baseline, 0.2*args.baseline], ylim, 'r--')
    AX[0].plot([1.2*args.baseline, 1.2*args.baseline], ylim, 'r--')
    AX[0].plot([x_crop[-1] - 0.2*args.baseline, x_crop[-1] - 0.2*args.baseline], ylim, 'r--')
    AX[0].plot([x_crop[-1] - 1.2*args.baseline, x_crop[-1] - 1.2*args.baseline], ylim, 'r--')
    AX[0].get_xaxis().set_visible(False)
    AX[0].get_yaxis().set_visible(False)
    AX[0].set_xlim([x_crop[0], x_crop[-1]])
    AX[0].set_ylim(ylim)

    # dynspec plot
    AX[1].imshow(ds_crop, aspect = 'auto', extent = [0, x_crop[-1], 0, 336])
    AX[1].plot([0.2*args.baseline, 0.2*args.baseline], [0,  336], 'r--')
    AX[1].plot([1.2*args.baseline, 1.2*args.baseline], [0, 336], 'r--')
    AX[1].plot([x_crop[-1] - 0.2*args.baseline, x_crop[-1] - 0.2*args.baseline], [0, 336], 'r--')
    AX[1].plot([x_crop[-1] - 1.2*args.baseline, x_crop[-1] - 1.2*args.baseline], [0, 336], 'r--')
    AX[1].set_xlabel("Time [ms]")
    AX[1].set_ylabel("Bandwidth [MHz]")

    fig.tight_layout()
    fig.subplots_adjust(hspace = 0)
    

    # save plot
    plt.savefig(f"{args.ofile.split('@')[0]}_bline_plot{label}.png")









def _proc(args, pol):
    """
    Main processing function
    """

    # initialise parameters
    sphase = None       # starting phase
    rbounds = None      # bounds for baseline correction
    chanflag = None     # channel flagging for baseline correction

    # flagging of channels should be done as soon as possible
    
    # loop through expected stokes suite
    stokes_params = "IQUV" if len(args.pols) > 1 else "I"

    for S in stokes_params:

        # make dynamic spectra
        ds = make_ds(pol['X'], pol['Y'], S, args.nFFT)

        if args.do_chanflag:
            ds_raw = ds.copy()

        # remove first channel (zero it), because reasons.... (probably due to instrumentation)
        ds[0] *= 1e-12

        # channel flagging
        if S == "I":
            chanflag, chanflag_known = flag_chan(ds, 10, 1000, args, rbounds)

        ## fold if a pulsar has been inputted
        if args.pulsar:
            ds, sphase = pulse_fold(ds, args.DM, args.cfreq, args.bw, args.MJD0, args.MJD1, 
                                      args.F0, args.F1, chanflag, sphase)
            
            if args.do_chanflag:
                ds_raw, sphase = pulse_fold(ds_raw, args.DM, args.cfreq, args.bw, args.MJD0, args.MJD1, 
                                      args.F0, args.F1, chanflag_known, sphase)
            
        if args.bline:

            ## get baseline corrections
            bs_mean, bs_std, rbounds = baseline_correction(ds, args.sigma, args.guard,
                                            args.baseline, args.tN, chanflag, rbounds)

            # check if failed
            if rbounds is None:
                return
        
            if args.do_chanflag:
                bs_mean_raw, bs_std_raw, rbounds = baseline_correction(ds_raw, args.sigma, args.guard,
                                             args.baseline, args.tN, chanflag_known, rbounds)
            
            ## Apply baseline corrections
            ds -= bs_mean[:, None]
            ds /= bs_std[:, None]

            if args.do_chanflag:
                ds_raw -= bs_mean_raw[:, None]
                ds_raw /= bs_std_raw[:, None]
                        
            # re-do channel flagging with proper baseline corrections, also plot bline diagnostics
            if S == "I":
                chanflag, _ = flag_chan(ds, 10, 1000, args, rbounds)
                plot_bline_diagnostic(ds, rbounds, chanflag, args)

                if args.do_chanflag:
                    plot_bline_diagnostic(ds_raw, rbounds, chanflag_known, args, "_nochanflag")

        
        if S == "I":
            np.save("flagged_channels.npy", chanflag)
            np.save("ASKAP_bad_channels.npy", chanflag_known)


        # flag channels in dynamic spectrum
        ds[chanflag] = np.nan

        ## save data
        print(f"Saving stokes {S} dynamic spectra...")
        np.save(args.ofile.replace("@", S), ds)

        if args.do_chanflag:
            ds_raw[chanflag_known] = np.nan
            print(f"Saving non-flagged stokes {S} dynamic spectra...")
            np.save(args.ofile.replace("@", f"{S}_nochanflag"), ds_raw)

        




if __name__ == "__main__":
    # main block of code

    ## get args
    args = get_args()


    ## load data
    pol = load_data(args)


    ## make dynamic spectra
    _proc(args, pol)


    print("Completed!")
