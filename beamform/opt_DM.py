from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from dedisperse import dedisperse
import numpy as np
import matplotlib.pyplot as plt

def _main():
    args = get_args()

    X = np.load(args.x)
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

    parser.add_argument("-x", type=str, help="X complex time series")
    parser.add_argument("-y", type=str, help="Y complex time series")
    parser.add_argument("--DM0", type=float, help="Baseline DM")
    parser.add_argument(
        "-d", "--minDM", type=float, default=0, help="DM range start"
    )
    parser.add_argument(
        "-D", "--maxDM", type=float, default=10, help="DM range end"
    )
    parser.add_argument(
        "-s", "--DMstep", type=float, default=0.01, help="DM range step"
    )
    parser.add_argument(
        "--dt", 
        type=int, 
        default=50, 
        help="Time resolution to average to in us"
    )
    parser.add_argument(
        "--f0",
        type=float,
        help="Central frequency in MHz"
    )
    parser.add_argument(
        "--bw",
        type=float,
        default=336,
        help="Bandwidth in MHz"
    )

    return parser.parse_args()


def running_mean(x, N):
    # https://stackoverflow.com/a/27681394
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def do_DM(X, Y, DM, dt, f0, bw):
    # do fft here because dedispersion operates in-place
    X_f = np.fft.fft(X)
    Y_f = np.fft.fft(Y)

    X_f_dd = dedisperse(X_f, DM, f0, bw)
    Y_f_dd = dedisperse(Y_f, DM, f0, bw)

    X_dd = np.fft.ifft(X_f_dd)
    Y_dd = np.fft.ifft(Y_f_dd)

    I = np.abs(X_dd)**2 + np.abs(Y_dd)**2
    I_red = running_mean(I, 336*dt)

    return np.max(I_red)


if __name__ == "__main__":
    _main()