# generates deripple coeffs

import numpy as np
from scipy import fftpack as fft
from scipy import io
from scipy.interpolate import interp1d


def generate_deripple(nfft, res, dir):
    N = 1536
    OS_De = 27.0
    OS_Nu = 32.0
    dr = dir + "/ADE_R6_OSFIR.mat"
    h = io.loadmat(dr)["c"][0]
    passbandLength = int(((nfft / 2) * OS_De) / OS_Nu)
    multiple = int(1536 / res)

    # Modifications: pre-pad h with zeros
    #                use scipy's fft instead of numpy's, as np's was segfaulting for huuuuge multiple*passbandLength*2
    h_0 = np.zeros(multiple * passbandLength * 2)
    h_0[: h.shape[0]] = h
    temp = abs(fft.fft(h_0))  # ,multiple*passbandLength*2))
    print(
        "saving {}".format(
            dir + "/deripple_res" + str(res) + "_nfft" + str(nfft)
        )
    )
    np.save(dir + "/deripple_res" + str(res) + "_nfft" + str(nfft), temp)


def _main():
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

    parser = ArgumentParser(
        description="Script description",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-n", "--nfft", type=int, help="FFT length", default=1000000
    )
    parser.add_argument("-r", "--res", type=int, default=6)
    values = parser.parse_args()

    try:
        print("FFT length: " + str(values.nfft))
        print("Resolution: " + str(values.res))
        generate_deripple(values.nfft, values.res)
    finally:
        print(
            "done generating deripple factors. Use interp1d("
            + str(values.res)
            + "*np.arange(len({OUTPUT})),{OUTPUT})"
        )


if __name__ == "__main__":
    _main()
