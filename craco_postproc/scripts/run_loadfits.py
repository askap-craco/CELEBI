#!/usr/bin/env python3

###############################################################################
# craftpy2 stage 2.3: loadfits
# Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
# Last modified: 2021-05-17
#
# Runs loadfits.py
###############################################################################

import os
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

prefix = "[run_loadfits.py]"

parser = ArgumentParser(
    description="Launch loadfits.py flexibly",
    formatter_class=ArgumentDefaultsHelpFormatter,
)
parser.add_argument("-i", "--info", help="Info file", required=True)
parser.add_argument(
    "-m",
    "--mode",
    help="Mode. Must be one of finder, field, " "or cal.",
    required=True,
)
args = parser.parse_args()

if args.mode not in ["finder", "field", "cal"]:
    print(f"{prefix} ERROR: mode {args.mode} invalid!")
    sys.exit()

# Parse info file into dict
info = {}
for line in open(args.info):
    spl = line[:-1].split("=")
    info[spl[0]] = spl[1]

cmd = "loadfits.py"
cmd += f" -u {os.getpid()}"  # use pid for user number to guarantee uniqueness
cmd += f' --antlist={info["antlist"]}'
cmd += " -s 27"

if args.mode == "finder":
    # Need to iterate over bins
    # Find number of bins from binconfig
    num_bins = -1
    for line in open(f"{args.mode}/craftfrb.finder.binconfig"):
        if "PULSAR" in line:
            spl = line[:-1].split(" ")
            num_bins = int(spl[-1])

    if num_bins < 0:
        print(f"{prefix} ERROR: Invalid number of bins ({num_bins})!")
        sys.exit()

    for bin in range(num_bins):
        bin_cmd = cmd
        bin_cmd += f' -f {args.mode}/{info["label"]}_{args.mode}_B{bin}.FITS'
        bin_cmd += f" -o {args.mode}_{bin}"
        bin_cmd += (
            f' {args.mode}/{info["label"]}_{args.mode}_CARD?_B{bin}.FITS'
        )

        print(bin_cmd)

else:
    # No need to iterate over bins
    cmd += f' -f {args.mode}/{info["label"]}_{args.mode}.FITS'
    cmd += f" -o {args.mode}"
    cmd += f' {args.mode}/{info["label"]}_{args.mode}_CARD?.FITS'

    print(f"{prefix} {cmd}")
    os.system(cmd)
