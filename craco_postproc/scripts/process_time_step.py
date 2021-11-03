#!/usr/bin/env python3

###############################################################################
# craftpy2 stage 2.2: process time step
# Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
# Last modified: 2021-05-19
#
# Runs processTimeStep.py
###############################################################################

import os
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from astropy.coordinates import SkyCoord as sc

prefix = "[process_time_step.py]"

parser = ArgumentParser(
    description="Launch processTimeStep.py flexibly",
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

cards = (1, 2, 3, 4, 5, 6, 7)

cmd = "processTimeStep.py"
cmd += f' -f {info["fcm"]}'
cmd += " -b 4"
cmd += " --card CARD"
cmd += " -k"
cmd += " --slurm"
cmd += f' --name={info["label"]}_{args.mode}'
cmd += f" -o data_{args.mode}"

if args.mode == "finder":
    # operating on FRB data
    cmd += f' -t {info["data_frb"]}'
    cmd += f' --ra {info["ra_frb"]} -d{info["dec_frb"]}'
    cmd += f' -i {info["int_time"]}'
    cmd += f" -p {args.mode}/craftfrb.finder.binconfig"

elif args.mode == "field":
    # operating on FRB data centred on the beam centre
    # get RA and Dec into hms and dms units
    pos_beam = sc(info["ra_beam"], info["dec_beam"], unit="deg")
    ra_beam_hms = pos_beam.ra.to_string("hour")
    dec_beam_dms = pos_beam.dec.to_string("deg")

    cmd += f' -t {info["data_frb"]}'
    cmd += f" --ra {ra_beam_hms} -d{dec_beam_dms}"

elif args.mode == "cal":
    # operating on calibrator data
    # get RA and Dec into hms and dms units
    cmd += f' -t {info["data_cal"]}'
    cmd += f' --ra {info["ra_cal"]} -d{info["dec_cal"]}'

cmd += f' > {info["label"]}/{args.mode}_cardCARD.log &'

for card in cards:
    card_cmd = cmd.replace("CARD", str(card))
    print(f"{prefix} {card_cmd}")
    os.system(card_cmd)
