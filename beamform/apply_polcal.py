##################################################
# Author:   Tyson Dial                           #
# Email:    tdial@swin.edu.au                    #
# Date (created):     20/11/2023                 #
# Date (updated):     18/12/2023                 #
##################################################
#  Apply POL CAL [apply_polcal.py]               #
#                                                #
# This script applyies derived pol cal solutions #
# to X and Y pol HTR FRB products                #
#                                                #
##################################################

# imports
import numpy as np
from scipy.fft import fft, ifft, next_fast_len
import argparse
from os import path


# empty class to store attributes
class empty_class():
    pass




def get_args():
    """
    Get arguments passed during script call

    """

    parser = argparse.ArgumentParser(
        description = "Apply Pol cal solutions"
    )

    ## data arguments
    parser.add_argument("-x", help = "X polarisation Time series", type = str)
    parser.add_argument("-y", help = "Y polarisation Time series", type = str)
    parser.add_argument("--soln", help = ".npy file containing solutions", type = str)


    ## observation arguments
    parser.add_argument("--cfreq", help = "Central frequency [MHz]", type = float)
    parser.add_argument("--bw", help = "Bandwidth [MHz]", type = float)


    ## output arguments
    # parser.add_argument("--odir", help = "Output directory for calibrated X,Y polarisations", type = str)
    # parser.add_argument("--ofile", help = "Prefix for output calibrated X pol files", type = str)
    parser.add_argument("--xout", help = "Output filename for X calib data", type = str)
    parser.add_argument("--yout", help = "Output filename for Y calib data", type = str)


    ## additional arguments
    parser.add_argument("--fast", help = "Apply zero-padding to fasten fft", action = "store_true")

    args = parser.parse_args()


    return args











def load_data(args):
    """
    Load in data

    #== inputs ==#
    args:           apply_polcal arguments

    #== outputs ==#
    X:              X polarisation time series
    Y:              Y polarisation time series
    par:            Pol cal solutions

    """
    
    # load in X, Y pol 
    X = np.load(args.x, mmap_mode = 'r')
    Y = np.load(args.y, mmap_mode = 'r')

    # load in polcal solutions
    par = empty_class()

    with open(args.soln, 'r') as f:
        par.psi = float(f.readline().split(':')[1])        # rotation offsert
        par.rm = float(f.readline().split(':')[1])         # rotation measure
        par.tau = float(f.readline().split(':')[1])        # time delay
        par.phi = float(f.readline().split(':')[1])        # phase delay
        par.Lfrac = float(f.readline().split(':')[1])      # Linear fraction free parameter
        par.Vfrac = float(f.readline().split(':')[1])      # Circular fraction free parameter
        par.alpha = float(f.readline().split(':')[1])      # ellipticity angle
        par.ellipticity = int(f.readline().split(':')[1])   # boolean, if ellipticity was modelled or not
    

    return X, Y, par










def apply_soln(args, X, Y, par):
    """
    Apply solutions

    #== inputs ==#
    args:           apply_polcal arguments
    X:              X polarisation
    Y:              Y polarisation
    par:            polcal solutions

    #== outputs ==#
    X_cal:              X polarisation (corrected)
    Y_cal:              Y polarisation (corrected)

    """
    print("[APPLY POLCAL]: Applying Polarisation leakage solutions...")

    ## Apply pol leakage solutions via Fourier Transform
    # FFT
    next_size = X.size
    if args.fast:
        next_size = next_fast_len(next_size)
        print(f"Zero-padding to speed up FFT: {X.size} -> {next_size}")
        print("Note: the more relative zero padding, the less accurate the results...")
    
    X_cal = fft(X, next_size)

    # apply rotation
    f_fine = np.linspace(args.cfreq - args.bw/2, args.cfreq + args.bw/2*(next_size/X.size), 
                         next_size)[::-1]
    X_cal *= np.exp(-2j * np.pi * f_fine * par.tau - par.phi * 1j)

    # iFFT
    X_cal = ifft(X_cal, next_size)[:X.size]

    # free up memory
    f_fine = None


    

    Y_cal = Y.copy()

    if par.ellipticity:

        print("[APPLY POLCAL]: Applying ellipticity correction...")

        # apply ellipticity correction
        Y_cal = -X_cal * 1j * np.sin(par.alpha) + Y * np.cos(par.alpha)
        X_cal = X_cal * np.cos(par.alpha) - Y * 1j * np.sin(par.alpha)

        Y = Y_cal.copy()


    print("[APPLY POLCAL]: Applying Polarisation rotation offset...")

    ## Apply rotational offset
    Y_cal = X_cal * np.sin(par.psi) + Y * np.cos(par.psi)
    X_cal = X_cal * np.cos(par.psi) - Y * np.sin(par.psi)


    return X_cal, Y_cal












# Entry point
if __name__ == "__main__":

    # get arguments
    args = get_args()

    # load data
    X, Y, par = load_data(args)

    # apply solutions
    X_cal, Y_cal = apply_soln(args, X, Y, par)

    # save data
    print("[APPLY POLCAL]: Saving calibrated X,Y data...")
    # np.save(path.join(args.odir, f"{args.ofile}_X.npy"), X_cal)       # save X polarisation
    # np.save(path.join(args.odir, f"{args.ofile}_Y.npy"), Y_cal)       # save Y polarisation
    np.save(args.xout, X_cal)   # save X polarisation
    np.save(args.yout, Y_cal)   # save Y polarisation


    print("[APPLY POLCAL]: Completed.")
