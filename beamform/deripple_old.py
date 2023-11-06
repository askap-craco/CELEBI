import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

import numpy as np
from scipy.interpolate import interp1d
from joblib import Parallel, delayed


def _main():
    start = time.time()
    args = get_args()
    sum_f = np.load(args.f)
    sum_f_deripple = deripple(sum_f, args.coeffs, args.l, args.bw, args.cpus)
    np.save(args.o, sum_f_deripple)
    end = time.time()
    print(f"deripple.py finished in {end-start} s")


def get_args() -> Namespace:
    """Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description="Deripples fine spectrum",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", help="Fine spectrum file")
    parser.add_argument("-l", type=int, help="FFT length")
    parser.add_argument("-o", help="Output file")
    parser.add_argument(
        "-c", "--coeffs", help="Deripple coefficients file"
    )
    parser.add_argument(
        "--bw", type=float, help="Spectrum bandwidth in MHz", default=336
    )
    parser.add_argument(
        "--cpus", type=int, help="Number of cpus to parallelise across", default=1
    )
    return parser.parse_args()


def deripple(
    FFFF: np.ndarray, coeffs_fname: str, fftLength: int, bw: float, cpus: int
) -> np.ndarray:
    """Deripple a fine spectrum.

    :param FFFF: Fine spectrum to be derippled
    :type FFFF: :class:`np.ndarray`
    :param coeffs_fname: Derippling coefficients filename
    :type coeffs_fname: str
    :param fftLength: Probably can be inferred from shape of FFFF
    :type fftLength: int
    :param bw: Bandwidth of spectrum in MHz
    :type bw: float
    :param cpus: Number of CPUs to parallelise across
    :type cpus: int
    :return: Derippled fine spectrum
    :rtype: :class:`np.ndarray`
    """
    print("derippling....")
    FFFF = FFFF[0, :, 0]

    # ASKAP Parameters
    N = 1536
    OS_De = 27.0
    OS_Nu = 32.0
    passbandLength = int(((fftLength / 2) * OS_De) / OS_Nu)

    coeffs = np.load(coeffs_fname, mmap_mode="r")

    print("Interpolating...")
    interp = interp1d(6 * np.arange(len(coeffs)), coeffs)

    print("Calculating deripple...")
    deripple = np.ones(passbandLength + 1) / abs(
        interp(np.arange(passbandLength + 1))
    )

    def do_deripple(chan, ii):
        FFFF[ii + chan * passbandLength * 2] = (
            FFFF[ii + chan * passbandLength * 2]
            * deripple[passbandLength - ii]
        )
        FFFF[passbandLength + ii + chan * passbandLength * 2] = (
            FFFF[passbandLength + ii + chan * passbandLength * 2]
            * deripple[ii]
        )        

    if cpus > 1:
        Parallel(n_jobs=cpus, require="sharedmem")(
            delayed(do_deripple)(chan, ii)
            for ii in range(passbandLength) for chan in range(bw)
        )
    else: 
        for chan in range(bw):
            for ii in range(passbandLength):
                do_deripple(chan, ii)

    return FFFF


if __name__ == "__main__":
    _main()
