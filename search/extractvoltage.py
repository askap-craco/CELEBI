##################################################
# Author: Apurba Bera                            #
##################################################
# Extract voltage data around the pulses         #          
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



def get_args():
    parser = ArgumentParser(description="Extract voltages around the pulse")

    # X voltage stream (Optional individually, but mutually required with -y)
    parser.add_argument("-x", type=str, default=None, help="X polarisation data file (optional, but at least one of -x or -y must be provided)")

    # Y voltage stream (Optional individually, but mutually required with -x)
    parser.add_argument("-y", type=str, default=None, help="Y polarisation data file (optional, but at least one of -x or -y must be provided)")

    # Location file
    parser.add_argument("--locfile", type=str, help="Location file")

    # half length of extracted data in us
    parser.add_argument("--hlen", type=float, help="half length in us")

    # Bandwidth in MHz
    parser.add_argument("--bwmhz", type=float, help="bandwidth in MHz")

    args = parser.parse_args()

    # Enforce that at least one polarisation is provided
    if args.x is None and args.y is None:
        parser.error(
            "\n\n"
            "CRITICAL ERROR: No polarisation data provided.\n"
            "You must provide at least one of the -x (X polarisation) or -y (Y polarisation) arguments.\n"
            " -> For single-polarisation data, provide just the one available.\n"
            " -> For dual-polarisation data, provide both."
        )

    return args




def load_data(args):
    
    # load XY data as memory map for memory efficiency. 
    # Conditionally load to support single pol.
    xx = np.load(args.x+".npy", mmap_mode = "r") if args.x else None
    yy = np.load(args.y+".npy", mmap_mode = "r") if args.y else None

    return xx,yy





def extract_htr(args, xx, yy):
    """
    Extract voltages around the pulse. 
    Dynamically supports single-polarisation or dual-polarisation data.
    """

    locdat  =   np.loadtxt(args.locfile)
    if( locdat.ndim < 2 ):
        locdat = np.array([locdat])

    locus   =   locdat[:,0]
    tresus  =   1.0/args.bwmhz

    # Determine reference array to check length bounds safely based on what is available
    ref_arr = xx if xx is not None else yy
    
    # This acts as a secondary safety net, though argparse handles the primary validation.
    if ref_arr is None:
        print("Error: No polarisation data (X or Y) could be loaded. Aborting extraction.")
        return 1

    for lo in range(0, len(locus)):
        lsamp   =   max(int((locus[lo] - args.hlen)/tresus), 0)
        rsamp   =   min(int((locus[lo] + args.hlen)/tresus), len(ref_arr))

        if xx is not None:
            xcrop   =   xx[lsamp:rsamp]
            print(f"X crop shape: {xcrop.shape}")
            np.save("{}_{}.npy".format("zoom_x",lo),xcrop)
        
        if yy is not None:
            ycrop   =   yy[lsamp:rsamp]
            print(f"Y crop shape: {ycrop.shape}")
            np.save("{}_{}.npy".format("zoom_y",lo),ycrop)

    return 0
    


if __name__ == "__main__":
    # main block of code

    # get arguments
    args = get_args()

    # get voltage data
    xx, yy = load_data(args)

    # Extract voltages
    extract_htr(args, xx, yy)
