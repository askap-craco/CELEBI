##################################################
# Author: Apurba Bera                            #
# Contributions: Based on work done by Danica    #
# Scott and Tyson Dial                           #
#                                                #
##################################################
# Search for additional pulses                   #          
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
    parser.add_argument("-s", type=str, required=True, help="Dynspecfile")

    # for correct frequency ranges
    parser.add_argument("-f", type=float, help="Central frequency (MHz)")

    # for labeling purposes
    parser.add_argument("-l", "--label", type=str, help="FRB label")

    # Time averaging
    parser.add_argument("--t_avgs", help = "Time resolutions to search",
                        nargs='+', default = [1, 10, 100, 1000, 10000], type = int)
    
    parser.add_argument("--thresh", help = "S/N threshold ",
                        default = 10.0, type = float)

    parser.add_argument("--fsub", help = "Number of frequency sub-bands",
                        default = 2, type = int)

    parser.add_argument("--suboverlap", help = "Overlap between subbands",
                        default = 50, type = int)      

    parser.add_argument("--edgech", help = "Edge channels to exclude",
                        default = 1, type = int)               

    return parser.parse_args()






def t_average(x, N):
    """
    Average x array by some factor N
    """
    if (N == 1):
        return x.copy()

    len_new = (x.shape[1] // N)

    return np.nanmean(x[:,:len_new*N].reshape(x.shape[0], len_new, N), axis = 2)






def load_data(args):
    
    # load dynspec as memory map for memory efficiency
    idspec = np.load(args.s, mmap_mode = "r")

    return idspec





def search_htr(args, ids):
    """
    Search additional bursts at different t resolutions 

    Parameters
    ----------
    args        - > function arguments
    ids         - > stokes I dynspec
    """

    t_arr = args.t_avgs
    nch   = ids.shape[0]
    ch0   = np.zeros(args.fsub, dtype=int)
    ch1   = np.ones(args.fsub, dtype=int) * nch
    wch   = (nch // args.fsub)

    # Configure sub-bands
    for ns in range(0, args.fsub):
        ch0[ns]  = (ns * wch) - (args.suboverlap // 2)
        ch1[ns]  = ch0[ns] + wch + (args.suboverlap // 2)

    ch0  = np.clip(ch0, args.edgech, nch - args.edgech)
    ch1  = np.clip(ch1, args.edgech, nch - args.edgech)

    candmask = np.zeros(ids.shape[1], dtype=int)

    # Configure sub-band time sries
    for ns in range(0, args.fsub):
        tI = np.nanmean(ids[ ch0[ns] : ch1[ns] ], axis = 0)

        # Time average
        for ta in t_arr:
            tI00  = t_average(tI.reshape(1, tI.size), ta).flatten()
            tI0   = tI00.copy()

            # Find robust RMS
            trms  = np.nanstd(tI0)
            tI0[tI00 > (args.thresh * trms)] = np.nan
            trms  = np.nanstd(tI0)
            tI0[tI00 > (args.thresh * trms)] = np.nan
            trms  = np.nanstd(tI0)
            dets  = np.nonzero( tI00 > (args.thresh * trms) )[0]

            print("Working on ta, tI00_shape, dets_shape: ", ta, tI00.shape, dets.shape)

            if(dets.shape[0] > 0):                
                print(dets)
                for det in dets:
                    candmask[ (det * ta) : (det+1) * ta ] = 1

            # Find periods
            fpow  = (np.abs(np.fft.rfft(tI00)))**2
            ffrq  = np.fft.rfftfreq(len(tI00), d=ta)
            np.savetxt("subband_"+str(ns)+"_"+"nt_"+str(ta)+".txt", np.column_stack((1.0e-3/ffrq,fpow)), fmt="%.3e  %.1e")

            fig   = plt.figure(figsize=(20.0, 8.0))
            plt.plot(1.0e-3/ffrq[5:], fpow[5:], 'b-', lw=1.0)
            plt.title("Channels "+str(ch0[ns])+" - "+str(ch1[ns])+" nt = "+str(ta))
            plt.xscale('log')
            #plt.yscale('log')
            plt.xlim([0.09, 1100.0])
            plt.xlabel("Period (ms)")
            plt.ylabel("Power")
            plt.savefig("subband_"+str(ns)+"_"+"nt_"+str(ta)+".png")
            plt.close()

    candedg   = np.ediff1d(candmask)
    candlr    = np.nonzero(candedg)[0] + 1

    cands     = np.cumsum(candmask)
    cands2    = np.ediff1d(cands[candlr])

    candl     = candlr[0::2]
    candw     = cands2[0::2]

    cands     = np.column_stack((candl,candw))
    np.savetxt("candfile.txt", cands, fmt="%ld  %ld")

    return cands
    




def plot_cands(args, ids, cands):
    """
    Plot candidate profiles and dynamic spectra
    """

    if (cands.ndim < 2):
        cands = cands.reshape(1,cands.size)

    for ci in range(0,len(cands)):

        cds     = ids[:,cands[ci,0] - 15*cands[ci,1]: cands[ci,0] + 17*cands[ci,1]]
        avds    = t_average(cds, 1 + (cands[ci,1] // 4))
        
        # Average 4 channels
        avds    = np.nanmean(np.reshape(avds, (84, 4, avds.shape[1])), axis=1)
        avds[:(args.edgech // 4)] = np.nan
        avds[-(args.edgech // 4):] = np.nan

        bkgns   = 1.48*np.nanmedian(np.abs(avds - np.nanmedian(avds)))

        cts     = np.nanmean(avds, axis=0)

        fig 	= plt.figure(figsize=(8.0,5.0))	
        ax 		= fig.add_axes([0.15, 0.50, 0.83,0.40])
        ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
        ax.axhline(c='c',ls='--',lw=0.5)	
        ax.plot(cts, 'b-', lw=1.0)	
        ax.set_xticklabels([])
        ax.set_xlim([0,len(cts)])
        ax.set_ylabel(r'Normalized intensity')
        ax.set_title("Cand "+str(ci)+": "+str(cands[ci,0])+" us -- width = "+str(cands[ci,1])+" us")

        ax1 		= 	fig.add_axes([0.15, 0.10, 0.83,0.40])
        ax1.imshow(avds, origin='lower', interpolation='none', aspect='auto', cmap='magma', vmin=-bkgns, vmax=5*bkgns)	
        ax1.set_xlim([0,len(cts)])
        ax1.set_xlabel(r'Time in units of '+str(1 + (cands[ci,1] // 4))+' us + ['+str(cands[ci,0] - 15*cands[ci,1])+' us]')
        ax1.set_ylabel(r'Channel')

        plt.savefig("cand_"+str(ci)+".png",	transparent=False, format='png', dpi=100)
        plt.close()
    
    
    








if __name__ == "__main__":
    # main block of code

    # get arguments
    args = get_args()

    # get burst starting time and load in stokes dynspecs
    ids = load_data(args)

    # Search dynamic spectrum
    cands = search_htr(args, ids)

    plot_cands(args, ids, cands)
