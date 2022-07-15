import numpy as np
from scipy import fft, io
import sys

res = 6

def _main():
    nfft = int(sys.argv[1])
    dr = sys.argv[2]    # path to ADE_R6_OSFIR.mat
    N = 1536
    OS_De = 27.0
    OS_Nu = 32.0
    h = io.loadmat(dr)["c"][0]
    passbandLength = int(((nfft / 2) * OS_De) / OS_Nu)
    multiple = int(1536 / res)

    h_0 = np.zeros(multiple * passbandLength * 2)
    h_0[: h.shape[0]] = h
    temp = abs(fft.fft(h_0))
    np.save(f"deripple_res{res}_nfft{nfft}.npy", temp)


if __name__ == "__main__":
    _main()
