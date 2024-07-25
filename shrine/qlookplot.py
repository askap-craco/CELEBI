#!/usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import os
import numpy as np
import matplotlib.pyplot as plt

"""
Make quick look plots of structure maximized dynamic spectrum
"""

def _main():

    args 	= get_args()

    i_ds 	= np.load(args.I, mmap_mode='r').T
    
    smdm	= np.load(args.S)

    f0 		= args.f0                      # centre frequency (MHz)
    bw 		= args.bw                      # bandwidth (MHz)
    dt 		= args.dt                      # time resolution to return in us
    dt0 	= args.dt0                 	   # Input time resolution in us
    
    nchan 	= i_ds.shape[1]
    freqs 	= get_freqs(f0, bw, nchan)
    
    Ism		= dedisperse_ic(i_ds, smdm[0], freqs, dt0)
    
    tavg	= int(np.round(dt/dt0))
    tlen	= int(Ism.shape[0]/tavg)
    
    Iavg	= np.zeros((tlen,nchan), dtype=float)
    for tt in range(0,tlen):
    	for cc in range(0,nchan):
    		Iavg[tt,cc]	=	np.nanmean(Ism[tt*tavg:tt*tavg+tavg,cc])
    
    I_ts	=	np.nanmean(Iavg, axis=1)
    
    fig 	= plt.figure(figsize=(5.0,5.0))	
    ax 		= fig.add_axes([0.15, 0.58, 0.83,0.40])
    ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
    ax.axhline(c='c',ls='--',lw=0.5)	
    ax.plot(I_ts, 'b-', lw=1.0, label="$\delta$DM = %.6f"%smdm[0])	
    ax.set_xticklabels([])
    ax.set_xlim([0,tlen])
    plt.legend(loc='upper right')
    ax.set_ylabel(r'Normalized flux density')

    ax1 		= 	fig.add_axes([0.15, 0.10, 0.83,0.48])
    ax1.imshow(Iavg.T, origin='lower', interpolation='none', aspect='auto', cmap='plasma')	
    ax1.set_xlim([0,tlen])
    ax1.set_xlabel(r'Time (%d us)'%dt)
    ax1.set_ylabel(r'Channel')

    plt.savefig("{}_idsmdm.png".format(args.label),	transparent=False, format='png', dpi=100)
    plt.close()

def get_args() -> ArgumentParser:
    """
    Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description = "Make quick look plot of structure maximized intensity profile", 
        formatter_class = ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-l", "--label",
        type = str,
        help = "FRB label"
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
        "-I",
        type = str,
        help = "Input intensity dynamic spectrum"
    )
    parser.add_argument(
        "-S",
        type = str,
        help = "SMDM result file"
    )
    return parser.parse_args()


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


if __name__ == "__main__":
    _main()
