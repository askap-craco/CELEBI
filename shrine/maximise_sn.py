from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct, idct #https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.dct.html

from dm_processing import get_kc

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

    delta_I=I_data-I_smooth # De-trended noise
    

    # S/N calculation

    # Calculation arguments
    window_step = 1
    min_window_size = max(window_step,1)
    max_window_size = int(10000/args.dt)

    # Obtain parameters for our noise
    # noise seems uniform across DMs, this is probably okay (?)
    constant_noise = np.min(I_smooth)
    noise_STD = np.std(delta_I)
    zeroed_I = I_smooth-constant_noise
    # Pre-calculate all possible denominators to save time
    time_roots = noise_STD*np.sqrt(np.arange(max_window_size))

    # Prep lists
    max_SNR_list = []
    max_SN_window_starts = []
    max_SN_window_lengths = []
    
    for iDM in zeroed_I:
        running_max_SN = -1
        for window_start in range(0,len(iDM)-min_window_size,window_step):
            S = np.cumsum(iDM[window_start:window_start+max_window_size])
            SN = S[min_window_size:]/time_roots[min_window_size:len(S)]
            max_SN = np.max(SN)
            if max_SN > running_max_SN:
                running_max_SN = max_SN
                max_SN_window_start = window_start
                max_SN_window_length = np.argmax(SN)
                
        max_SNR_list.append(running_max_SN)
        max_SN_window_starts.append(max_SN_window_start)
        max_SN_window_lengths.append(max_SN_window_length)
    
    max_SN_index = np.argmax(max_SNR_list)

    plot_SN_v_DM(DM_data, max_SNR_list, args)
    plot_max_sn_window(I_smooth, max_SN_index, max_SN_window_starts, max_SN_window_lengths, args)

    plot_I_at_max(I_data[max_SN_index], I_smooth[max_SN_index], args)

    # Find uncertainty bounds
    max_SN = np.max(max_SNR_list)
    upp_bound_idx = -1
    low_bound_idx = -1

    i = max_SN_index
    while i < len(max_SNR_list):
        i += 1
        if max_SN - max_SNR_list[i] >= 1:
            break
    upp_bound_idx = i

    i = max_SN_index
    while i >= 0:
        i -= 1
        if max_SN - max_SNR_list[i] >= 1:
            break
    low_bound_idx = i


    # S/N of DENOISED signal
    denoised_signal_sum = np.sum(zeroed_I[max_SN_index][max_SN_window_starts[max_SN_index]:max_SN_window_starts[max_SN_index]+max_SN_window_lengths[max_SN_index]])
    denoised_signal_noise = np.concatenate((zeroed_I[:max_SN_window_starts[max_SN_index]],zeroed_I[max_SN_window_lengths[max_SN_index]:]))
    denoised_signal_SN = denoised_signal_sum/(np.std(denoised_signal_noise)*np.sqrt(max_SN_window_lengths[max_SN_index]))

    #print(MaxSNRs)
    summary_file = open(f"{args.label}_SN_summaryfile.txt", "w")
    summary_file.write(f"//begin maximise_sn summary//\n/*\n")
    summary_file.write(f"Max SNR at index {max_SN_index}\n")
    summary_file.write(f"Uncertainty bounds at indexes {low_bound_idx} - {upp_bound_idx}\n")
    summary_file.write(f"Corresponds to a delta DM of {DM_data[max_SN_index]} +{DM_data[upp_bound_idx]-DM_data[max_SN_index]} -{DM_data[max_SN_index]-DM_data[low_bound_idx]}\n")
    summary_file.write(f"Calculated with window from {max_SN_window_starts[max_SN_index]} to {max_SN_window_starts[max_SN_index]+max_SN_window_lengths[max_SN_index]}\n")
    summary_file.write(f"Corresponds to a time window from {max_SN_window_starts[max_SN_index]*args.dt/1000}us to {(max_SN_window_starts[max_SN_index]+max_SN_window_lengths[max_SN_index])*args.dt/1000}us\n")
    summary_file.write(f"Max SNR was {np.max(max_SNR_list)}\n")
    summary_file.write(f"*/\n//end maximise_sn summary//\n\n\n")
    summary_file.close()

    print(f"S/N of the ALREADY DENOISED signal was {denoised_signal_SN}")

    SN_delta = DM_data - DM_data[max_SN_index]
    np.savetxt(f"{args.label}_SN_DMs.dat", SN_delta)

    if args.save_all:
        plot_SN_v_DM(DM_data, max_SNR_list, args)

        plot_max_sn_window(I_smooth, max_SN_index, max_SN_window_starts, max_SN_window_lengths, args)

        plot_I_at_max(I_data[max_SN_index], I_smooth[max_SN_index], args)

        np.save(f"{args.label}_I_smooth.npy", I_smooth)
        DM_file = open(f"DM.txt", "w")
        DM_file.write(f"{max_SN_index}")
        DM_file.close()


def get_args() -> ArgumentParser:
    """
    Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description="Calculates signal to noise ratio maximising delta DM for pre-dedispersed data.", 
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

def plot_SN_v_DM(DM_data: np.ndarray, max_SNR_list: np.ndarray, args):
    plt.plot(DM_data, max_SNR_list)
    plt.grid(color='k', linestyle='--', linewidth=0.5)
    plt.xlabel(r"$\Delta$DM ($\mathregular{pc\ cm^{-3}}$)")
    plt.ylabel('Signal to Noise Ratio')
    plt.savefig(f"{args.label}_SN.png")
    plt.clf()

def plot_max_sn_window(I_smooth: np.ndarray, max_SN_index: int, max_SN_window_starts: np.ndarray, max_SN_window_lengths: np.ndarray, args):
    plt.plot([i*args.dt/1000 for i in range(len(I_smooth[max_SN_index]))], I_smooth[max_SN_index])
    plt.xlabel("Time (ms)")
    plt.ylabel("I")
    plt.axvspan(
        max_SN_window_starts[max_SN_index]*args.dt/1000,
        (max_SN_window_starts[max_SN_index]+max_SN_window_lengths[max_SN_index])*args.dt/1000,
        color = 'red',
        alpha = 0.3
        )
    plt.savefig(f"{args.label}_Max_SN_window.png")
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
    plt.plot([i*args.dt/1000 for i in range(len(I_series))], I_series, alpha=0.5, color='red')
    plt.plot([i*args.dt/1000 for i in range(len(I_series))], I_smooth_series, color='red')
    plt.xlabel("Time (ms)")
    plt.ylabel("I")
    plt.savefig(f"{args.label}_{args.dt}us_I_SN_max.png")
    plt.clf() 

if __name__ == "__main__":
    _main()
