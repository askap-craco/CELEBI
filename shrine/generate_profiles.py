#!/usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
#from tqdm import tqdm

"""
Originally written by Timothy Perrett
Modified by Danica Scott to get profiles via dynamic spectra and allow
variable crop duration
Further modified by AB to allow 
    masking of bad channels (26June2024)
    directly passing HTR dynamic spectra
"""

def _main():

    args = get_args()

    Iact = np.load(args.I, mmap_mode='r')
    I    = Iact.T    
    
    f0 = args.f0                      # centre frequency (MHz)
    bw = args.bw                      # bandwidth (MHz)
    dDM = args.dDM                    # DM step
    cDM = args.cDM                    # DM count
    dt = args.dt                      # time resolution to return in us
    dt0 = args.dt0                    # Input time resolution in us
    if cDM == 0:
        DMs = np.arange(args.DM_low, args.DM_high, dDM)  # DM range
    else:
        DMs = np.linspace(args.DM_low, args.DM_high, cDM)
    
    I_profs = []
    
    nchan = I.shape[1]
    
    print("Found %d channels and %d samples"%(nchan,I.shape[0]))
    
    freqs = get_freqs(f0, bw, nchan)

    for DM in DMs:
        I_profs.append(do_DM(I, DM, freqs, dt0, dt, bw))

    I_profs = np.array(I_profs)
    #I_profs = np.transpose(np.transpose(I_profs)[3000:8000])
        
    np.save(f"{args.label}_DMs.npy", DMs)
    np.save(f"{args.label}_I_{dt}us.npy", I_profs)

    np.savetxt(f"{args.label}_DMs.dat", DMs)
    np.savetxt(f"{args.label}_I_{dt}us.dat", I_profs)

    summary_file = open(f"{args.label}_profile_summaryfile.txt", "w")
    summary_file.write(f"//begin generate_profiles summary//\n/*\n")
    #summary_file.write(f"X and Y files were of length {len(X)}\n")
    summary_file.write(f"DM range was from {DMs[0]} to {DMs[-1]}\n")
    summary_file.write(f"DM range was made up of {len(DMs)} steps of size {DMs[1]-DMs[0]}\n")
    summary_file.write(f"I profile was of shape {I_profs.shape}\n")
    summary_file.write(f"*/\n//end generate_profiles summary//\n\n\n")
    summary_file.close()

def get_args() -> ArgumentParser:
    """
    Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description = "Generates an I(DM,t) profile for a given FRB.", 
        formatter_class = ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-l", "--label",
        type = str,
        help = "FRB label"
    )
    parser.add_argument(
        "-d", "--DM",
        type = str,
        help = "Dispersion measure to which the input has been dedispersed in pc/cm3"
    )
    parser.add_argument(
        "-r", "--DM_range",
        type = float,
        help = "Range of delta DMs either side of 0 to profile in pc/cm3",
        default = 5
    )
    parser.add_argument(
        "-L", "--DM_low",
        type = float,
        help = "Lower bound of delta DM range",
        default = -5
    )
    parser.add_argument(
        "-H", "--DM_high",
        type = float,
        help = "Upper bound of delta DM range",
        default = 5
    )
    parser.add_argument(
        "--dDM",
        type = float,
        help = "Delta DM step size in pc/cm3",
        default = 0.1
    )
    parser.add_argument(
        "--cDM",
        type = int,
        help = "Number of DM steps to divide the range into",
        default = 0
    )
    parser.add_argument(
        "-t0", "--dt0",
        type = int,
        help = "Input time resolution to return in us",
        default = 1
    )
    parser.add_argument(
        "-t", "--dt",
        type = int,
        help = "Time resolution to return in us",
        default = 1
    )
    parser.add_argument(
        "-f", "--f0", 
        type = float,
        help = "Central frequency of observation in MHz"
    )
    parser.add_argument(
        "-b", "--bw",
        type = int,
        help = "Bandwidth of observation in MHz",
        default = 336
    )
    parser.add_argument(
        "--crop_dur",
        type = int,
        help = "Duration of the crop in ms",
        default=100
    )
    parser.add_argument(
        "-I",
        type = str,
        help = "I profile filename (for doing the crop)"
    )
    parser.add_argument(
        "--force_peak",
        type = int,
        help = "Force the peak to be at this time step (in ms)",
        default=None
    )
    return parser.parse_args()


def do_DM(I0_ds, DM, freqs, dt0, dt, bw):    

    I_ds = dedisperse_ic(I0_ds, DM, freqs, dt0)

    #	AB: Replacing mean by median and avoiding normalization
    for i in range(I_ds.shape[1]):
        I_ds[:,i] -= np.nanmedian(I_ds[:,i])
	
	#	Identify and flag bad channels
    chanmads	=	np.zeros(I_ds.shape[1], dtype=float)
	
    for i in range(I_ds.shape[1]):
	    chanmads[i]	=	np.nanmedian(np.abs(I_ds[:,i]))
	
    medmads		=	np.nanmedian(chanmads)
    devmads		=	chanmads - medmads
    madmads		=	np.nanmedian(np.abs(devmads))
    badcount	=	0
	
    for i in range(I_ds.shape[1]):
        if(devmads[i] > 6*madmads):
            I_ds[:,i]	=	np.nan
            badcount	+=	1
    
    print("Saving intensity dynamic spectrum (%d x %d) in %s \n"%(I_ds.shape[0],I_ds.shape[1],os.getcwd()))
    np.save("dspec_i_"+str(DM)+".npy",I_ds)
    
    print("Flagged %d channels...\n"%badcount)
    print("Excluding 5+5 edge channels\n")
	
    #	Add channels to consturct time profiles
    I_dd 	= np.nansum(I_ds[:,5:-5],axis=1)
	
    ndt 	= int(np.round(dt/dt0)) 
    print("Tscrunching %d samples\n"%ndt)
    I_red 	= tscrunch(I_dd, ndt)
    np.save("dts_i_"+str(DM)+".npy",I_red)

    return I_red


def dedisperse_ic(I0_ds, DM, freqs, dtus):
    """
    Incoherent de-dispersion
    
    Inputs
        I0_ds   - Input Stokes I dynamic spectrum (order -- time, frequency)
        DM      - DM to correct for
        freqs   - Frequency array in MHz
        dtus    - Time resolution in us
        
    Output
        I_ds    - De-dispersed dynamic spectrum
    """
    k_dm    = 4.149
    f_ref   = np.nanmedian(freqs)/1.0e3         #   Converting to GHz
    fghz    = freqs/1.0e3

    I_ds    = np.zeros(I0_ds.shape, dtype=float)
    
    for cc in range(0,I0_ds.shape[1]):
        delshft     = - int(np.round(( k_dm * DM * ( (1.0/fghz[cc]**2) - (1.0/f_ref**2) ) ) *1.0e3 / dtus))
        I_ds[:,cc]  = np.roll(I0_ds[:,cc], delshft)
    
    return I_ds


def get_freqs(f0: float, bw: float, nchan: int) -> np.ndarray:
    """
    Create array of frequencies.

    The returned array is the central frequency of `nchan` channels
    centred on `f0` with a bandwidth of `bw`.

    :param f0: Central frequency (arb. units, must be same as `bw`)
    :type f0: float
    :param bw: Bandwidth (arb. units, must be same as `f0`)
    :type bw: float
    :param nchan: Number of channels
    :type nchan: int
    :return: Central frequencies of `nchan` channels centred on `f0`
        over a bandwidth `bw`
    :rtype: :class:`np.ndarray`
    """
    fmin = f0 - bw / 2
    fmax = f0 + bw / 2

    chan_width = bw / nchan

    freqs = np.linspace(fmax, fmin, nchan, endpoint=False) + chan_width / 2

    return freqs


def tscrunch(a: np.ndarray, n: int) -> np.ndarray:
    """
    Reduces the time resolution of array `a` by a factor `n`.

    :param a: The array to scrunch
    :type a: :class:`np.ndarray`
    :param n: The factor by which to reduce `a`
    :type n: int
    :return: Array with reduced time resolution
    :rtype: :class:`np.ndarray`
    """

    a_red = []

    for i in range(int(a.shape[0] / n)):
        a_red.append(np.nansum(a[i * n:(i + 1) * n], axis=0))

    a_red = np.array(a_red) / np.sqrt(n)

    return a_red


def generate_dynspec(
        t_ser: np.ndarray, nchan: int = 336,
    ) -> np.ndarray:
    """
    Creates a dynamic spectrum at the highest time resolution from the 
    given time series.

    :param t_ser: input time series of voltages
    :param nchan: number of frequency channels [Default = 336]
    :return: dynamic spectrum of voltages
    :rtype: :class:`np.ndarray`
    """
    dynspec = np.zeros(
        (int(t_ser.shape[0] / nchan), nchan), dtype=np.complex64
    )

    for i in range(int(t_ser.shape[0] / nchan)):
        dynspec[i, :] = np.fft.fft(t_ser[i * nchan : (i + 1) * nchan])

    return dynspec


if __name__ == "__main__":
    _main()
