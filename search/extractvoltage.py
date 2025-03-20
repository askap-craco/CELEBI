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

    # X voltage stream
    parser.add_argument("-x", type=str, required=True, help="X data")

    # Y voltage stream
    parser.add_argument("-y", type=str, required=True, help="Y data")

    # Location file
    parser.add_argument("--locfile", type=str, help="Location file")

    # half length of extracted data in us
    parser.add_argument("--hlen", type=float, help="half length in us")

    # Bandwidth in MHz
    parser.add_argument("--bwmhz", type=float, help="bandwidth in MHz")

    return parser.parse_args()




def load_data(args):
    
    # load XY data as memory map for memory efficiency
    xx = np.load(args.x+".npy", mmap_mode = "r")
    yy = np.load(args.y+".npy", mmap_mode = "r")

    return xx,yy





def extract_htr(args, xx, yy):
    """
    Extract voltages around the pulse
    """

    locdat  =   np.loadtxt(args.locfile)
    if( locdat.ndim < 2 ):
        locdat = np.array([locdat])

    locus   =   locdat[:,0]
    tresus  =   1.0/args.bwmhz

    for lo in range(0, len(locus)):
        lsamp   =   max(int((locus[lo] - args.hlen)/tresus), 0)
        rsamp   =   min(int((locus[lo] + args.hlen)/tresus), len(xx))
        xcrop   =   xx[lsamp:rsamp]
        ycrop   =   yy[lsamp:rsamp]

        print(xcrop.shape,ycrop.shape)

        np.save("{}_{}.npy".format("zoom_x",lo),xcrop)
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
