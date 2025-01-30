from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np
import matplotlib.pyplot as plt

def _main():

    args = get_args()

    DM_data = np.load(args.DM_data)
    I_data = np.load(args.I_data)
    smooth_data = np.load(args.smooth_data)

    plt.plot([i*args.dt/1000 for i in range(smooth_data.shape[1])], smooth_data[args.i_st], label=f"Structure Maximised\n($\Delta$ DM = {round(DM_data[args.i_st], 5)})", color='blue')
    plt.plot([i*args.dt/1000 for i in range(smooth_data.shape[1])], smooth_data[args.i_sn], '--', label=f"S/N Maximised\n($\Delta$ DM = {round(DM_data[args.i_sn], 5)})", color='red')
    plt.xlabel("Time (ms)")
    plt.ylabel("I")
    plt.legend()
    plt.savefig(f"{args.label}_structure_SN_comparison.png")
    plt.clf()

def get_args() -> ArgumentParser:
    """
    Parse command line arguments

    :return: Command line argument parameters
    :rtype: :class:`argparse.Namespace`
    """
    parser = ArgumentParser(
        description="Plots a comparison of the smoothed signal at the structure and S/N maximising delta DMs.", 
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
        "-i_st",
        type=int,
        help="Index of the structure maximising delta DM"
    )
    parser.add_argument(
        "-i_sn",
        type=int,
        help="Index of the signal to noise ratio maximising delta DM"
    )
    parser.add_argument(
        "-d", "--DM_data",
        type=str,
        help="Path to the DM data.",
    )
    parser.add_argument(
        "-I", "--I_data",
        type=str,
        help="Path to the I data.",
    )
    parser.add_argument(
        "-s", "--smooth_data",
        type=str,
        help="Path to the smoothed I data.",
    )
    return parser.parse_args()

if __name__ == "__main__":
    _main()
