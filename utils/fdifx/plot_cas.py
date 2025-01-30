# imports
from fdifx import difx
import numpy as np
import matplotlib.pyplot as plt 
from util import get_inputfiles, get_sort_index
import argparse




def get_args():
    """
    Get script args
    """

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-d", help = "Data difx directory (where all the c*_f* dirs are).", type = str)

    return parser.parse_args()






def plot_cas(args):
    """
    Plot CAS

    Parameters
    ----------
    args : argparse.argumentparser
        arguments for script
    """

    # get inputs files
    inputfiles = get_inputfiles(args.d)


    # get freqs
    freq = []
    for i, inp in enumerate(inputfiles):
        idifx = difx(inp)
        freq += [idifx.freq_info.freq[0]]
        freq += [idifx.freq_info.freq[1]]
    
    freq = np.array(freq)

    tbins = idifx.nbins
    fbins = len(inputfiles) * 2     # assume 4MHz coarse channels


    # array for CAS 
    cas = np.zeros((fbins, tbins), dtype = float)


    # loop over inputfiles and run get_cas
    for i, inp in enumerate(inputfiles):

        print(f"Progress [cas]: {(i+1)/(len(inputfiles)):.3%}")

        idifx = difx(inp)

        # loop over each bin
        for j in range(1, tbins):
            idifx.loadbin(bin = j)
            cas[2*i,j] = idifx.get_cas(idifx.freq_info.bandsActive[0])
            cas[2*i+1,j] = idifx.get_cas(idifx.freq_info.bandsActive[1])
    
    
    # sort cas frequencies
    fsort_idx = get_sort_index(freq)
    cas = cas[fsort_idx]


    # plot cas
    fig, AX = plt.subplots(1, 1, figsize = (14,10))
    tbins -= 1
    AX.imshow(cas[::-1,1:], aspect = 'auto', extent = [ 0, tbins, freq[0], freq[-1]])
    AX.set_xlabel("time bin", fontsize = 16)
    AX.set_ylabel("Freq [MHz]", fontsize = 16)

    # adjust figure
    fig.tight_layout()
    fig.subplots_adjust(hspace = 0)

    plt.savefig("cas.png")
    plt.show()






if __name__ == "__main__":
    # main code block


    # get args
    args = get_args()


    # plot CAS
    plot_cas(args)


    # complete
    print("Complete!")





