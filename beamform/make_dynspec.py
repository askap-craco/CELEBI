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
def average(x, axis = 0, N = 1):
    """
    scrunch data in either f or t

    ##==== inputs ====##
    x:              data
    N:              scrunch factor
    axis:           axis to scrunch in

    ##==== outputs ====## 
    x:              scrunched data

    """
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
    parser.add_argument("-x", help = "X pol time series", type = str)
    parser.add_argument("-y", help = "Y pol time series", type = str)
    parser.add_argument("--nFFT", help = "Number of frequency channels for final dynspec", 
                        type = int, default = 336)
    parser.add_argument("--bline", help = "Apply baseline correction", action = "store_true")


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

















def load_data(xfile, yfile):
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
    pol = {}

    pol['X'] = np.load(xfile, mmap_mode = 'r')
    pol['Y'] = np.load(yfile, mmap_mode = 'r')


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
    nsamps  = xpol.size                  # original number of samples in loaded dataset
    fnsamps = (nsamps // nFFT) * nFFT    # number of samples after chopping 
    nwind   = fnsamps // nFFT            # number of fft windows along time series

    # memeory block paramters
    nwinb = int(BLOCK_SIZE // (nFFT * BIT_SIZE))    # num windows in BLOCK
    nsinb = int(nwinb * nFFT)                       # num samples in BLOCK
    nblock= int(nwind // nwinb)                     # num BLOCKS

    # create empty array
    ds = np.zeros((nFFT, nwind), dtype = np.float32)

    b_arr = np.empty((0,2), dtype = int)
    for i in range(nblock):
        b_arr = np.append(b_arr, [[i*nwinb,(i+1)*nwinb]], axis = 0)
    # append extra block at end
    if nblock*nsinb < nsamps:
        b_arr = np.append(b_arr, [[(i+1)*nwinb,nwind]], axis = 0)

    # loop over blocks
    for i, b in enumerate(b_arr): # b is bounds of block in nFFT windows
        sb = b * nFFT
        wind_w = b[1] - b[0]
        ds[:,b[0]:b[1]] = Stk_Func[S](fft(xpol[sb[0]:sb[1]].copy().reshape(wind_w, nFFT), axis = 1),
                                 fft(ypol[sb[0]:sb[1]].copy().reshape(wind_w, nFFT), axis = 1)).T
        
        # print progress
        print(f"[MAKING DYNSPEC]:    [Progress] = {(i+1)/(nblock+1)*100:3.3f}%:    " + prog_str,
              end = '         \r')

    print("[MAKING DYNSPEC]:    [Progress] = 100.00%:    " + prog_str + "        \n")
    print(f"Made Dynamic spectra with shape [{ds.shape[0]}, {ds.shape[1]}]")

    return ds







def pulse_fold(ds, DM, cfreq, bw, MJD0, MJD1, F0, F1, sphase = None):
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

    Returns:
        ds_f (ndarray): Raw folded Dynamic spectra
        sphase (float): Starting phase of folding from original dynspec

    """
    print("Pulse Folding Dynspec...")
    # normalise needed to help find fold bounds
    ds_mean = np.mean(ds, axis = 1)[:, None]
    ds_std = np.std(ds, axis = 1)[:, None]

    ## Calculate Period T in [s]
    T = 1/(F0 - F1 * (MJD1 - MJD0)*86400)
    print(f"with period T = {T}")
    dt = 1e-6 * (ds.shape[0]/336) # get time resolution of dynspec

    ## Fold dynamic spectra
    fold_w = int(T / dt)          # fold width in samples (assumed dt = 1 us)
    fold_n_init = int(ds.shape[1]/fold_w)     # initial number of folds

    # get dispersion sweep, calculate number of "broken" pulse periods
    # due to dipsersion.
    top_band = args.cfreq + bw/2
    bot_band = args.cfreq - bw/2
    DM_sweep = 4.14938e3 * DM * (1/bot_band**2 - 1/top_band**2) # DM sweep in seconds
    P_sweep = int(DM_sweep/T) + 1
    print(f"DM sweep: {DM_sweep} [ms]")
    print(f"Culling {P_sweep} Periods to the left due to DM sweep")



    # find index of peak in second period, then get starting phase
    if sphase is None:
        search_crop = ds[:,P_sweep*fold_w:(P_sweep + 1)*fold_w]
        search_crop = (search_crop - ds_mean)/ds_std
        pulse2i = np.mean(search_crop, axis = 0).argmax()
        pulse2i += P_sweep * fold_w
        sphase = pulse2i - int(fold_w/2)
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















def baseline_correction(ds, sigma: float = 5.0, guard: float = 1.0, 
                        baseline: float = 50.0, tN: int = 50, rbounds = None):
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
    # rough normalisation needed to find burst
    ds_n = (ds - np.mean(ds, axis = 1)[:, None])/(np.std(ds, axis = 1)[:, None])

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
        ds_r = average(ds_n, axis = 1, N = tN)
        rmean = np.mean(ds_r, axis = 1)
        rstd = np.std(ds_r, axis = 1)

        ds_rn = ds_r - rmean[:, None]
        ds_rn /= rstd[:, None]

        
        ## find burst bounds
        print("Looking for bounds of burst...")
        # get peak, crop rms and do rough S/N calculation
        t_rn = np.mean(ds_rn, axis = 0)
        peak = np.argmax(t_rn)
        rms_w = get_units_avg(baseline)
        rms_crop = np.roll(t_rn, int(rmsg * ds_rn.shape[1]))[peak-rms_w:peak+rms_w]
        rms = np.mean(rms_crop**2)**0.5

        # calculate S/N
        t_sn = t_rn / rms
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






def plot_bline_diagnostic(ds, rbounds, args):
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

    # get time series
    t_crop = np.mean(ds_crop, axis = 0)

    # get time axis in ms/ or 1000 x dt
    x_crop = np.linspace(0, t_crop.size/1000*args.tN, t_crop.size)


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
    plt.savefig(f"{args.ofile.split('@')[0]}_bline_plot.png")









def _proc(args, pol):
    """
    Main processing function
    """

    # initialise parameters
    sphase = None       # starting phase
    rbounds = None      # bounds for baseline correction
    
    # loop over full stokes suite
    for S in "IQUV":

        # make dynamic spectra
        ds = make_ds(pol['X'], pol['Y'], S, args.nFFT)

        # remove first channel (zero it)
        ds[0] *= 1e-12

        ## fold if a pulsar has been inputted
        if args.pulsar:
            ds, sphase = pulse_fold(ds, args.DM, args.cfreq, args.bw, args.MJD0, args.MJD1, 
                                      args.F0, args.F1, sphase)

        if args.bline:
            ## get baseline corrections
            bs_mean, bs_std, rbounds = baseline_correction(ds, args.sigma, args.guard,
                                            args.baseline, args.tN, rbounds)

            if S == "I":
                plot_bline_diagnostic(ds, rbounds, args)


            ## Apply baseline corrections
            ds -= bs_mean[:, None]
            ds /= bs_std[:, None]

        ## save data
        print(f"Saving stokes {S} dynamic spectra...")
        # np.save(f"{args.ofile}_{S}.npy", ds)
        np.save(args.ofile.replace("@", S), ds)




if __name__ == "__main__":
    # main block of code

    ## get args
    args = get_args()


    ## load data
    pol = load_data(args.x ,args.y)


    ## make dynamic spectra
    _proc(args, pol)


    print("Completed!")




