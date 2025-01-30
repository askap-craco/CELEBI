import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

import numpy as np
from joblib import Parallel, delayed


########################################################################
########################################################################
# UPDATED 08/09/2023 - TYSON DIAL                                      #
# Email: tdial@swin.edu.au                                             #
#                                                                      #
#                                                                      #
# UPDATE NOTE:                                                         #
# implemented vectorization, using numpy to apply deripple coeffs.     #
# Also updated scipy.interp -> np.interp1                              #
# Removed Scipy module                                                 #
#                                                                      #
#                                                                      #
#                                                                      #
#                                                                      #
########################################################################
########################################################################


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
        "--bw", type=int, help="Spectrum bandwidth in MHz", default=336
    )
    parser.add_argument(
        "--cpus", type=int, help="Number of cpus to parallelise across", default=1
    )
    return parser.parse_args()


def deripple(
    FFFF: np.ndarray, coeffs_fname: str, fftLength: int, bw: int, cpus: int
) -> np.ndarray:
    """Deripple a fine spectrum.

    :param FFFF: Fine spectrum to be derippled
    :type FFFF: :class:`np.ndarray`
    :param coeffs_fname: Derippling coefficients filename
    :type coeffs_fname: str
    :param fftLength: Probably can be inferred from shape of FFFF
    :type fftLength: int
    :param bw: Bandwidth of spectrum in MHz
    :type bw: int
    :param cpus: Number of CPUs to parallelise across
    :type cpus: int
    :return: Derippled fine spectrum
    :rtype: :class:`np.ndarray`
    """

    # bandwidth, i.e. number of coarse channels
    bw = int(bw) # MAKE SURE IT'S AN INTEGER
    
    # FFFF is the stitched 336 (or some other bw) MHz fine channel complex beamformed data. 
    FFFF = FFFF[0, :, 0]    
    
    # load in coefficients
    coeffs = np.load(coeffs_fname, mmap_mode="r")

    # calculate the number of fine channels in each coarse channel, which should simply be the size of the
    # buffer in samples / bw in MHz
    nfine = FFFF.size // bw

    # passbandLength is nfine/2 samples, however, whether nfine is odd or even, will change how
    # we create the deripple coeff array.
    isdiv2 = False
    if nfine % 2 == 0:
        passbandLength = nfine // 2
        isdiv2 = True

    else:   # if not divisible by two, extend passbandlength by one and truncat mid point when concatenating
        passbandLength = (nfine // 2) + 1
    
    print("Interpolating...")
    interp_x = np.interp(np.arange(passbandLength),6*np.arange(coeffs.size),coeffs)


    print("Calculating deripple...")
    deripple = np.ones(passbandLength) / np.abs(interp_x)

    #PRINT SOME INFOMATION
    print("FFT LENGTH (oversamp): {:d}".format(fftLength))
    print("PASS BAND LENGTH: {:d}".format(passbandLength))
    print("BAND WIDTH (MHz): {:d}".format(bw))

    print("derippling....")
    
    if isdiv2: # can simply mirror [deripple] to make coeffs array
        deripple = np.concatenate((deripple[::-1],deripple), axis = 0)
    else:   # truncate end sample 
        deripple = np.concatenate((deripple[:0:-1], deripple), axis = 0)

    print(f"number of samples per coarse channel: {nfine}")
    print(f"interpolated deripple coeff sample number: {deripple.size}")

    #reshape
    FFFF = FFFF.reshape(bw, nfine)
    
    #apply deripple
    FFFF *= deripple

    print("derippling Done.")
    #redo reshapping, make single fine spectrum
    return FFFF.flatten()

if __name__ == "__main__":
    _main()
