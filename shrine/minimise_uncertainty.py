from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.fftpack import dct, idct #https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.dct.html
from scipy import linalg #https://docs.scipy.org/doc/scipy/reference/linalg.html

from dm_processing import get_kc, uncertainty_calc, get_ranges_above_max

def _main():

    args = get_args()

    DM_data = np.load(f"{args.label}_DMs.npy"); #print(DMdata)
    I_data = np.load(f"{args.label}_I_{args.dt}us.npy"); #print(Idata)

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
    kc_set = range(max(10,real_kc-100),real_kc+100)

    results = []
    uncertainty_lower = []
    uncertainty_upper = []

    for kc in kc_set:
        #filter response
        fL=1/(1+(k/kc)**(2*O))

        #Pass DCT data through the combined filter to calculate structure parameter
        fL_diag=np.diag(fL) #make low-pass Filter into a diagonal matrix
        LPF_data=fL_diag@np.transpose(CI_data) #pass data through LPF
        I_smooth=idct(np.transpose(LPF_data), norm='ortho') #smooth data

        #1st derivative "high-pass" filter
        pi=np.pi
        hp=np.sqrt(2-2*np.cos((k-1)*pi/k_length)) #square root of the eigenvalues of D_1^T D_1

        #combined response
        filter=hp*fL
        filter_diag=np.diag(filter) #make Filter into a diagonal matrix

        #pass DCT data through the combined filter
        CI_filtered=filter_diag@np.transpose(CI_data)
        norm_CI_filtered=linalg.norm(CI_filtered, axis=0) #calculate structrue parameter, note numpy axis (see https://www.sharpsightlabs.com/blog/numpy-axes-explained/)


        #Uncertainty calculations
        delta_I=I_data-I_smooth #de-trended noise

        #calculate uncertainty proper
        uncertainty = uncertainty_calc(delta_I, LPF_data, filter_diag)

        #Second round of uncertainty calculations
        max_index = np.argmax(norm_CI_filtered)
        delta_delta_I=delta_I-delta_I[max_index] 

        relative_uncertainty = uncertainty_calc(delta_delta_I, LPF_data, filter_diag)

        #Find uncertainty in structure maximizing DM
        max_structure_parameter = norm_CI_filtered[max_index]

        adjusted_SPs = norm_CI_filtered + (norm_CI_filtered * relative_uncertainty) 

        possible_max_ranges = get_ranges_above_max(max_structure_parameter, adjusted_SPs)
        
        if possible_max_ranges[0][0] == 0:
            uncertainty_lower.append(DM_data[0]-0.1) 
        else:
            uncertainty_lower.append(DM_data[possible_max_ranges[0][0]])
        if len(possible_max_ranges[-1])==2:
            uncertainty_upper.append(DM_data[possible_max_ranges[-1][1]])
        else:
            uncertainty_upper.append(DM_data[-1]+0.1)
        results.append(DM_data[max_index])

    if args.save_all:
        plt.rcParams['font.family'] = 'cm'
        plt.rcParams['mathtext.fontset'] = 'cm'
        plt.rcParams['font.size'] = 14
        plt.plot(kc_set, results, color='black', label='Structure Maximising\nDelta DM')
        plt.plot(kc_set, uncertainty_lower, '--', color='black', label='Uncertainty Bounds')
        plt.plot(kc_set, uncertainty_upper, '--', color='black')
        plt.ylim(DM_data[0], DM_data[-1])
        plt.xlabel(r"$k_c$")
        plt.ylabel(r"$\rm{\Delta\ DM}\ (\rm{pc\ cm}^{-3})$")
        #plt.xlim(right = 70)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{args.label}_kc_var_uncertainty.png")
        plt.clf()

    np.savetxt(f"{args.label}_uncertainly_kc.dat", np.vstack([kc_set, uncertainty_lower, results, uncertainty_upper]).T)

    uncertainty_ranges = [uncertainty_upper[i]-uncertainty_lower[i] for i in range(len(kc_set))]

    min_uncertainty = uncertainty_ranges[0]
    min_indices = []
    for i, uncertainty_range in enumerate(uncertainty_ranges):
        if uncertainty_range < min_uncertainty:
            min_indices = [i]
            min_uncertainty = uncertainty_range
        elif uncertainty_range == min_uncertainty:
            min_indices.append(i)
    
    min_index = min_indices[int(len(min_indices)/2)] # middle value of range
    kc_file = open(f"kc.txt", "w")
    kc_file.write(f"{kc_set[min_index]}")
    kc_file.close()

    summary_file = open(f"{args.label}_uncertainty_summaryfile.txt", "w")
    summary_file.write(f"//begin minimise_uncertainty summary//\n/*\n")
    if args.force_kc is None:
        summary_file.write(f"kc was calculated to be {real_kc}\n")
    else:
        summary_file.write(f"kc was forced to be {real_kc}\n")
    summary_file.write(f"kc range tested was from {kc_set[0]} to {kc_set[-1]}\n")
    summary_file.write(f"smallest uncertainty range found was of size {min_uncertainty}\n")
    summary_file.write(f"this range had DM bounds {uncertainty_lower[min_index]} to {uncertainty_upper[min_index]}\n")
    summary_file.write(f"the number of kc values which produced a range of this size was {len(min_indices)}\n")
    summary_file.write(f"the selected uncertainty minimising kc value was {kc_set[min_index]}\n")
    summary_file.write(f"*/\n//end maximise_structure summary//\n\n\n")
    summary_file.close()

def get_args() -> ArgumentParser:
    """
    Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description="Calculates the value of kc which minimises uncertainty. Also generates a plot of structure maximising delta DM and uncertainty against kc.", 
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
        "-s", "--save_all",
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
