import time

import numpy as np
try:
    from scipy.fft import ifft
except:
    from scipy.fftpack import ifft


def _main():
    start = time.time()
    args = get_args()
    f = load(args.f)
    t = do_ifft(f)
    save(t, args.o)
    end = time.time()
    print(f"ifft.py finished in {end-start} s")


def get_args():
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

    parser = ArgumentParser(
        description="Performs ifft on given spectrum to obtain time series"
    )
    parser.add_argument("-f", help="Spectrum file to ifft")
    parser.add_argument("-o", help="Output file to save time series to")
    return parser.parse_args()


def load(fname):
    print(f"Loading {fname}")
    return np.load(fname)


def do_ifft(f):
    print("IFFTing")
    return ifft(f)


def save(t, fname):
    print(f"Saving {fname}")
    np.save(fname, t)


if __name__ == "__main__":
    _main()
