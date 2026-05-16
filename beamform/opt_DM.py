import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from dedisperse import dedisperse
import numpy as np
import matplotlib.pyplot as plt

def _main():
    args = get_args()

    expected_pols = args.pols
    X = None
    Y = None

    # Check expectations and load
    if 'X' in expected_pols:
        if not args.x or not os.path.exists(args.x):
            raise FileNotFoundError("CRITICAL ERROR: Polarisation 'X' was expected based on --pols, but X data file is missing or not provided.")
        X = np.load(args.x)
        
    if 'Y' in expected_pols:
        if not args.y or not os.path.exists(args.y):
            raise FileNotFoundError("CRITICAL ERROR: Polarisation 'Y' was expected based on --pols, but Y data file is missing or not provided.")
        Y = np.load(args.y)

    DMs = np.arange(args.minDM, args.maxDM+args.DMstep, args.DMstep)

    peaks = []
    for DM in DMs:
        peaks.append(do_DM(X, Y, DM, args.dt, args.f0, args.bw))
    
    print(DMs[np.argmax(peaks)]+args.DM0)

    plt.plot(DMs+args.DM0, peaks)
    plt.axvline(DMs[np.argmax(peaks)]+args.DM0)
    plt.xlabel("DM (pc/cm3)")
    plt.ylabel("max(I)")
    plt.tight_layout()
    plt.savefig("opt_DM.png")


def get_args():
    parser = ArgumentParser(
        "Optimise DM for S/N", 
        formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-x", type=str, default=None, help="X complex time series")
    parser.add_argument("-y", type=str, default=None, help="Y complex time series")
    parser.add_argument("--pols", nargs='+', required=True, help="List of expected polarisations (e.g., X Y or just X)")
    
    parser.add_argument("--DM0", type=float, help="Baseline DM")
    parser.add_argument("-d", "--minDM", type=float, default=0, help="DM range start")
    parser.add_argument("-D", "--maxDM", type=float, default=10, help="DM range end")
    parser.add_argument("-s", "--DMstep", type=float, default=0.01, help="DM range step")
    parser.add_argument("--dt", type=int, default=50, help="Time resolution to average to in us")
    parser.add_argument("--f0", type=float, help="Central frequency in MHz")
    parser.add_argument("--bw", type=float, default=336, help="Bandwidth in MHz")

    return parser.parse_args()


def running_mean(x, N):
    # https://stackoverflow.com/a/27681394
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def do_DM(X, Y, DM, dt, f0, bw):
    # Process X if available
    if X is not None:
        X_f = np.fft.fft(X)
        X_f_dd = dedisperse(X_f, DM, f0, bw)
        X_dd = np.fft.ifft(X_f_dd)
        I = np.abs(X_dd)**2
        
    # Process Y if available
    if Y is not None:
        Y_f = np.fft.fft(Y)
        Y_f_dd = dedisperse(Y_f, DM, f0, bw)
        Y_dd = np.fft.ifft(Y_f_dd)
        if X is not None:
            I += np.abs(Y_dd)**2
        else:
            I = np.abs(Y_dd)**2

    I_red = running_mean(I, int(bw*dt))
    return np.max(I_red)


if __name__ == "__main__":
    _main()
