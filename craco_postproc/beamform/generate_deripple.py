import numpy as np
from scipy import fft, io
import sys

def _main():
    nfft = int(sys.argv[1])
    dr = sys.argv[2]    # path to ADE_R6_OSFIR.mat
    N = 1536
    res = 6
    OS_De = 27.0
    OS_Nu = 32.0
    h = io.loadmat(dr)["c"][0]
    passbandLength = int(((nfft / 2) * OS_De) / OS_Nu)
    multiple = int(N / res)

    h_0 = np.zeros(multiple * passbandLength * 2)
    h_0[: h.shape[0]] = h
    np.save("h_0.npy", h_0)
    del h_0
    h_0 = np.load("h_0.npy", mmap_mode="r")

    # Save memory by computing rfft, then casting it into shape of fft
    # output
    temp = rfftout_to_fftout(fft.rfft(h_0, overwrite_x=True))
    np.save("temp", temp)
    del temp
    temp = np.load("temp.npy", mmap_mode="r")
    abstemp = np.abs(temp)
    np.save(f"deripple_res{res}_nfft{nfft}.npy", abstemp)


def rfftout_to_fftout(x):
    """
    Cast the output of an rfft into the shape it would be if you had done an fft
    """
    n = x.shape[0]
    fftout = np.zeros(2*n-2, dtype=complex)
    fftout[:n] = x
    fftout[n:] = np.conjugate(x[1:-1][::-1])
    return fftout



if __name__ == "__main__":
    _main()
