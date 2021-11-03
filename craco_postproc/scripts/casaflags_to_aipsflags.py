#!/usr/bin/env python3

###############################################################################
# casaflags_to_aipsflags
# Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
# Last modified: 2021-05-26
#
# Convert flags as determined in casa into the format required by AIPS
###############################################################################

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np

parser = ArgumentParser(
    description="Converts flags as determined in casa into "
    "the format required by AIPS",
    formatter_class=ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-c",
    "--casa",
    help="Comma-separated file containing flags"
    " as determined in casa. Each line "
    "should contain antenna number, start "
    "channel, and end channel. Set antenna"
    " to -1 to apply flag to all antennas",
)
parser.add_argument(
    "-o", "--outfile", help="File to write AIPS-formatted flags" " to."
)
args = parser.parse_args()

flags = np.loadtxt(args.casa, delimiter=",", dtype="int")

with open(args.outfile, "w") as f:
    f.write("dtimrang = 1  timeoff = 0\n\n")

    basestr = (
        "antennas=ANT bchan=START echan=END timerang=0,0,0,0,0,23,59,59 "
        "reason='RFI'/\n"
    )

    for flag in flags:
        # CASA is 0-indexed, AIPS is 1-indexed, so we need to increase ant and
        # start by 1. End is increased by 2 because AIPS doesn't include the
        # end channel. Also, in AIPS antennas=0 flags all antennas
        ant = flag[0] + 1
        start = flag[1] + 1
        end = flag[2] + 2

        flag_str = basestr.replace("ANT", str(ant))
        flag_str = flag_str.replace("START", str(start))
        flag_str = flag_str.replace("END", str(end))

        f.write(flag_str)
