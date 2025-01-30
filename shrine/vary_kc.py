from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.fftpack import dct, idct #https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.dct.html
from scipy import linalg #https://docs.scipy.org/doc/scipy/reference/linalg.html

from dm_processing import get_kc


def _main():

    args = get_args()

    DM_data = np.load(f"{args.label}_DMs.npy")
    I_data = np.load(f"{args.label}_I_{args.dt}us.npy")

    #Setup complete, start doing the math
    CI_data=dct(I_data, norm='ortho') #note "norm=ortho" to match MATLAB's dct
    dm_length,k_length=CI_data.shape

    #Low pass filter
    if args.force_kc == None:
        real_kc = get_kc(CI_data) #cutoff k index
    else:
        real_kc = args.force_kc
    O=3; #filter order
    k=np.linspace(1,k_length,k_length)
    kc_set = [int(real_kc/16),int(real_kc/8),int(real_kc/4),int(real_kc/2),real_kc,real_kc*2,real_kc*4,real_kc*6,real_kc*8]
    results = []
    data_sets = []
    for kc in kc_set:
        #filter response
        fL=1/(1+(k/kc)**(2*O))

        #Pass DCT data through the combined filter to calculate structure parameter
        fL_diag=np.diag(fL) #make low-pass Filter into a diagonal matrix
        LPF_I_spec=fL_diag@np.transpose(CI_data) #pass data through LPF
        I_smooth=idct(np.transpose(LPF_I_spec), norm='ortho') #smooth data

        #1st derivative "high-pass" filter
        pi=np.pi
        hp=np.sqrt(2-2*np.cos((k-1)*pi/k_length)) #square root of the eigenvalues of D_1^T D_1

        #combined response
        filter=hp*fL
        filter_diag=np.diag(filter) #make Filter into a diagonal matrix

        #pass DCT data through the combined filter
        CI_filtered=filter_diag@np.transpose(CI_data)
    
        norm_CI_filtered=linalg.norm(CI_filtered, axis=0) #calculate structrue parameter, note numpy axis (see https://www.sharpsightlabs.com/blog/numpy-axes-explained/)

        max_DM_index = np.argmax(norm_CI_filtered)
        
        data_sets.append(norm_CI_filtered)
        results.append(DM_data[max_DM_index])

    #single plot
    plt.style.use('dark_background')
    cmap = mpl.cm.get_cmap('bwr')
    kc_norm = mpl.colors.TwoSlopeNorm(vmin=np.min(kc_set),vmax=np.max(kc_set),vcenter=real_kc)
    lines = []
    for i, data in enumerate(data_sets):
        lines.append(plt.plot(DM_data,data,'--',color=cmap(kc_norm(kc_set[i])),label=kc_set[i]))
        plt.vlines(DM_data[np.argmax(data)],0,max(data),color=cmap(kc_norm(kc_set[i])),zorder=-i)
    lines[kc_set.index(real_kc)][0].set(linestyle='-')

    plt.legend(title = 'kc')
    plt.ylim(bottom=0)
    plt.xlabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.ylabel("Structure Parameter")
    
    plt.savefig(f"{args.label}_varied_kc_SPs.png")


def get_args() -> ArgumentParser:
    """
    Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description="Generates a plot showing the effects of varying kc.", 
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-l", "--label",
        type=str,
        help="FRB label"
    )
    parser.add_argument(
        "-t", "--dt",
        type=str,
        help="Time resolution of inputs in us.",
        default=1
    )
    parser.add_argument(
        "-d", "--DM",
        type=float,
        help="Dispersion measure to which the input has been dedispersed in pc/cm3"
    )
    parser.add_argument(
        "-s", "--saveAll",
        action="store_true",
        default=False,
        help="Save all plots.",
    )
    parser.add_argument(
        "-kc", "--force_kc",
        default=None,
        type=int,
        help="Value of kc to use.",
    )
    return parser.parse_args()
    

if __name__ == "__main__":
    _main()
