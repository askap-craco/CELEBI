#!/usr/bin/env python

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from astropy import units as u
from astropy.coordinates import SkyCoord as sc

parser = ArgumentParser(
    description="Execute jmfit on a given fits image at "
    + "the provided coordinates."
)
parser.add_argument("--askap", "-a", help="file of ASKAP positions")
parser.add_argument("--nvss", "-n", help="file to input NVSS positions to")
args = parser.parse_args()

print(
    "For each source the ASKAP file prints, copy/paste the RA and Dec\n"
    + "exactly as they appear into the NVSS catalog search. Find the source\n"
    + "in the list given by NVSS, and copy/paste its entry (BOTH LINES)\n"
    "into the prompt here."
)

NVSS_file = open(args.nvss, "w")
basestr = "{0},{1},{2},{3}\n"

for line in open(args.askap):
    RA, RA_err, Dec, Dec_err = tuple(line[:-1].split(","))
    print("RA:  {}".format(RA.replace(":", " ")))
    print("Dec: {}".format(Dec.replace(":", " ")))

    NVSS_vals = eval(input("> "))
    NVSS_errs = eval(input())

    NVSS_RA = NVSS_vals[:11].replace(" ", ":").replace("::", ":0")
    NVSS_Dec = NVSS_vals[12:23].replace(" ", ":").replace("::", ":0")
    NVSS_RA_err = float(NVSS_errs[7:11])
    NVSS_Dec_err = float(NVSS_errs[20:23])

    NVSS_RA_err_as = (
        sc(f"00:00:{NVSS_RA_err}", 0, unit="hour,deg").ra.to(u.arcsec).value
    )

    NVSS_file.write(
        basestr.format(NVSS_RA, NVSS_RA_err_as, NVSS_Dec, NVSS_Dec_err)
    )

NVSS_file.close()

"""
example NVSS catalog item
22 39  0.48 -16 00 57.7  0.56  130.6  16.0 <14.4 -43.8     18.62 -25.4 C2240M16  569.21  509.11
       0.03         0.6    11    4.7   2.8         2.7      0.46   0.5
"""
