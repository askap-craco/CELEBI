from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct, idct #https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.dct.html
from scipy import linalg #https://docs.scipy.org/doc/scipy/reference/linalg.html

from dm_processing import get_kc, uncertainty_calc, get_ranges_above_max


def _main():

    args = get_args()

    DM_data = np.load(f"{args.label}_DMs.npy")
    I_data = np.load(f"{args.label}_I_{args.dt}us.npy")

    #Setup complete, start doing the math
    CI_data=dct(I_data, norm='ortho') #note "norm=ortho" to match MATLAB's dct
    dm_length,k_length=CI_data.shape
    
    #Low pass filter
    if args.force_kc is None:
        kc = get_kc(CI_data) #cutoff k index
    else:
        kc = args.force_kc
    O=3 #filter order
    k=np.linspace(1,k_length,k_length)
    #filter response
    fL=1/(1+(k/kc)**(2*O))

    #Pass DCT data through the combined filter to calculate structure parameter
    fL_diag=np.diag(fL) #make low-pass Filter into a diagonal matrix
    LPF_data=fL_diag@np.transpose(CI_data) #pass data through LPF
    I_smooth=idct(np.transpose(LPF_data), norm='ortho') #smooth data

    np.save(f"{args.label}_I_smooth_{args.dt}us.npy", I_smooth)

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
    delta_DM_at_max = DM_data[max_index]

    adjusted_SPs = norm_CI_filtered + (norm_CI_filtered * relative_uncertainty) 

    possible_max_ranges = get_ranges_above_max(max_structure_parameter, adjusted_SPs)

    if len(possible_max_ranges) >= 1:
        if len(possible_max_ranges) == 1:
            one_range = True
        else:
            one_range = False
        # possible_max_ranges ALWAYS finds a minimum, so this is safe
        min_delta_DM = DM_data[possible_max_ranges[0][0]]
        if len(possible_max_ranges[-1]) == 2:
            max_delta_DM = DM_data[possible_max_ranges[-1][1]]
        else:
            max_delta_DM = None
    else:
        # Something has gone wrong if this runs
        # possible_max_ranges collects everything where `SP` >= `SP at max`
        # so at a minimum the max range should consist of `SP_at_max`
        min_delta_DM = None
        max_delta_DM = None
        one_range = None

    np.savetxt(f"{args.label}_SPs.dat", norm_CI_filtered)
    np.savetxt(f"{args.label}_Uncertainties.dat", uncertainty)
    np.savetxt(f"{args.label}_Relative_Uncertainties.dat", relative_uncertainty)

    if args.save_all:

        DM_file = open(f"DM.txt", "w")
        DM_file.write(f"{max_index}")
        DM_file.close()

        plot_DM_index(DM_data, args)

        # Intensity-time plots
        plot_I_at_max(I_data[max_index],I_smooth[max_index],args)
        plot_noise_at_max(delta_I[max_index], args)

        # I(DM,t) plots
        plot_noisy_I_DM_t(I_data, DM_data, args)
        plot_smooth_I_DM_t(I_smooth, DM_data, args)
        plot_detrended_noise(delta_I, DM_data, args)
        plot_relative_detrended_noise(delta_delta_I, DM_data, args)
        plot_DCT_spectrum(CI_data, args, kc)

        # Uncertainty-time plots
        plot_uncertainty(uncertainty, DM_data, args)
        plot_relative_uncertainty(relative_uncertainty, DM_data, args)

        # Structure parameter plots
        plot_SP(norm_CI_filtered, DM_data, args)
        plot_adjusted_SP(adjusted_SPs, DM_data, args)
	
    np.save(f"{args.label}_smdm_result.npy",np.array([delta_DM_at_max, nonetoflzero(min_delta_DM)-delta_DM_at_max, nonetoflzero(max_delta_DM)-delta_DM_at_max])) 

    #Write a nice tidy file
    summary_file = open(f"{args.label}_structure_summaryfile.txt", "w")
    summary_file.write(f"//begin maximise_structure summary//\n/*\n")
    summary_file.write(f"FRB Label: {args.label}\t//should match params.label\n")
    summary_file.write(f"Time Resolution: {args.dt}us\t//should match params.timescale\n")
    summary_file.write(f"Saving: {args.save_all}\t//should match params.saving\n")
    summary_file.write(f"kc: {kc}\n")
    if args.force_kc is None:
        summary_file.write(f"\tForced kc: {False}\n")
    else:
        summary_file.write(f"\tForced kc: {True}\n\t(This could be because of a provided force_kc value or because minimise_uncertainty was used.\n")
    summary_file.write(f"Initial DM: {args.DM}\n")
    summary_file.write(f"Structure Maximising Delta DM: {delta_DM_at_max}\n")
    summary_file.write(f"Uncertainty in Structure Maximising Delta DM: ")
    if min_delta_DM is not None:
        if min_delta_DM != np.min(DM_data):
            if max_delta_DM is not None:
                summary_file.write(f"{min_delta_DM-delta_DM_at_max}/+{max_delta_DM-delta_DM_at_max}\n")
            else:
                summary_file.write(f"{min_delta_DM-delta_DM_at_max}/+unknown\n")
                summary_file.write(f"\tUpper bound on uncertainty exceeds delta DM range.\n")
        else:
            if max_delta_DM is not None:
                summary_file.write(f"{min_delta_DM-delta_DM_at_max}/+{max_delta_DM-delta_DM_at_max}\n")
                summary_file.write(f"\tLower bound on uncertainty is equal to lower bound on delta DM range.\n")
                summary_file.write(f"\tThis probably means it is actually lower than {np.min(DM_data)}.\n")
            else:
                summary_file.write(f"{min_delta_DM-delta_DM_at_max}/+unknown\n")
                summary_file.write(f"\tLower bound on uncertainty is equal to lower bound on delta DM range.\n")
                summary_file.write(f"\tThis probably means it is actually lower than {np.min(DM_data)}.\n")
                summary_file.write(f"\tUpper bound on uncertainty exceeds delta DM range.\n")
    else:
        summary_file.write(f"-unknown/+unknown\n")
        summary_file.write(f"\tNo uncertainty range found. Something has gone wrong!\n")
    if not one_range:
        summary_file.write(f"\tUncertainty range was not continuous.\n")
        summary_file.write(f"\tRanges found were from:\n")
        for possible_range in possible_max_ranges:
            if len(possible_range) == 2:
                #good range
                summary_file.write(f"\t{DM_data[possible_range[0]]} to {DM_data[possible_range[1]]}\n")
            else:
                #bad range
                summary_file.write(f"\t{DM_data[possible_range[0]]} onwards. (Range finishes out of bounds).\n")
    summary_file.write(f"*/\n//end maximise_structure summary//\n\n\n")
    summary_file.close()


def get_args() -> ArgumentParser:
    """
    Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description="Calculates structure maximising delta DM for pre-dedispersed data.", 
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-l", "--label",
        type=str,
        help="FRB label"
    )
    parser.add_argument(
        "-t", "--dt",
        type=int,
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


def nonetoflzero(val):
    return float(0 if val is None else val)


def plot_DM_index(dm_data: np.ndarray, args):
    """
    Generates and saves a plot of delta DM vs the DM index.
    
    Will be a straight line unless someone goes out of their way to generate a non-linear delta DM range.
    
    :param DMdata: The set of DMs corresponding to the rows of Idata.
    :type DMdata: :class:`np.ndarray`
    :param args: Command line arguments.
    :type args: :class:`argparse.Namespace`
    """
    plt.plot(dm_data)
    plt.xlabel('DM Index')
    plt.ylabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.savefig(f"{args.label}_{args.dt}us_DM_index.png")
    plt.clf() 


def plot_noisy_I_DM_t(I_data: np.ndarray, DM_data: np.ndarray, args):
    """
    Generates and saves a plot of I(DM,t).
    
    :param Idata: The I(DM,t) profile generated by generate_profiles.py.
    :type Idata: :class:`np.ndarray`
    :param DMdata: The set of DMs corresponding to the rows of `Idata`.
    :type DMdata: :class:`np.ndarray`
    :param args: Command line arguments.
    :type args: :class:`argparse.Namespace`
    """
    imgplot_Inoisy = plt.imshow(I_data, aspect="auto", extent=(0, I_data.shape[1]*args.dt/1000, DM_data[-1], DM_data[0]))
    plt.xlabel("Time (ms)")
    plt.ylabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.savefig(f"{args.label}_{args.dt}us_I_DM_t.png")
    plt.clf() 

def plot_I_at_max(I_series: np.ndarray, I_smooth_series: np.ndarray, args):
    """
    Generates and saves a plot of the intensity time series at the structure maximising DM.
    
    :param I_series: The I(DM,t) profile generated by generate_profiles.py.
    :type I_series: :class:`np.ndarray`
    :param I_series: The smoothed I(DM,t) profile.
    :type I_series: :class:`np.ndarray`
    :param args: Command line arguments.
    :type args: :class:`argparse.Namespace`
    """
    plt.plot([i*args.dt/1000 for i in range(len(I_series))], I_series, alpha=0.5, color='blue')
    plt.plot([i*args.dt/1000 for i in range(len(I_series))], I_smooth_series, color='blue')
    plt.xlabel("Time (ms)")
    plt.ylabel("I")
    plt.savefig(f"{args.label}_{args.dt}us_I_max.png")
    plt.clf() 

def plot_DCT_spectrum(CI_data: np.ndarray, args, kc):
    plt.plot(np.transpose(np.abs(CI_data)),'.')
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel(r'$\mathregular{k}$ Index')
    plt.ylabel(r'DCT coefficients, $\mathregular{|C^Ti|}$')
    plt.axvline(x=kc)
    plt.grid(color='k', linestyle='--', linewidth=0.5)
    plt.savefig(f"{args.label}_{args.dt}us_DCT_v_k_index.png")
    plt.clf() 


def plot_smooth_I_DM_t(I_smooth: np.ndarray, DM_data: np.ndarray, args):
    imgplot_Ismooth = plt.imshow(I_smooth, aspect="auto", extent=(0, I_smooth.shape[1]*args.dt/1000, DM_data[-1], DM_data[0]))
    plt.xlabel("Time (ms)")
    plt.ylabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.savefig(f"{args.label}_{args.dt}us_I_DM_t_smoothed.png")
    plt.clf() 


def plot_SP(SP_data: np.ndarray, DM_data: np.ndarray, args):
    plt.plot(DM_data,SP_data,'-')
    plt.xlabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.ylabel('Structure Parameter')
    plt.grid(color='k', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(f"{args.label}_{args.dt}us_structure_parameter.png")
    plt.clf() 


def plot_detrended_noise(delta_I: np.ndarray,DM_data: np.ndarray, args):
    imgplot_DeltaI = plt.imshow(delta_I, aspect="auto", extent=(0, delta_I.shape[1]*args.dt/1000, DM_data[-1], DM_data[0]))
    plt.xlabel("Time (ms)")
    plt.ylabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.savefig(f"{args.label}_{args.dt}us_detrended_noise.png")
    plt.clf() 


def plot_relative_detrended_noise(delta_delta_I: np.ndarray, dm_data: np.ndarray, args):
    imgplot_DeltaDeltaI = plt.imshow(delta_delta_I, aspect="auto", extent=(0, delta_delta_I.shape[1]*args.dt/1000, dm_data[-1], dm_data[0]))
    plt.xlabel("Time (ms)")
    plt.ylabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.savefig(f"{args.label}_{args.dt}us_relative_detrended_noise.png")
    plt.clf() 


def plot_noise_at_max(noise_at_max: np.ndarray, args): 
    plt.plot([time/1000 for time in range(len(noise_at_max*args.dt))], noise_at_max)
    plt.xlabel("Time (ms)")
    plt.ylabel("I")
    plt.savefig(f"{args.label}_{args.dt}us_noise_max.png")
    plt.clf() 


def plot_uncertainty(uncertainty: np.ndarray, DM_data: np.ndarray, args):
    plt.plot(DM_data,100*uncertainty,'-')
    plt.xlabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.ylabel('Uncertainty (%)')
    plt.grid(color='k', linestyle='--', linewidth=0.5)
    plt.savefig(f"{args.label}_{args.dt}us_uncertainty.png")
    plt.clf() 


def plot_relative_uncertainty(relative_uncertainty: np.ndarray, DM_data: np.ndarray, args):
    plt.plot(DM_data,100*relative_uncertainty,'-')
    plt.xlabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.ylabel('Relative Uncertainty (%)')
    plt.grid(color='k', linestyle='--', linewidth=0.5)
    plt.savefig(f"{args.label}_{args.dt}us_relative_uncertainty.png")
    plt.clf()


def plot_adjusted_SP(adjusted_SPs: np.ndarray, DM_data: np.ndarray, args):
    plt.plot(DM_data,adjusted_SPs,'-')
    plt.xlabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.ylabel('Relative Uncertainty Adjusted Structure Parameter')
    plt.grid(color='k', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(f"{args.label}_{args.dt}us_adjusted_SP.png")
    plt.clf()

if __name__ == "__main__":
    _main()
