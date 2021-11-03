#!/usr/bin/env python

import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

parser = ArgumentParser(
    description="Execute jmfit on a given fits image at "
    + "the provided coordinates."
)
parser.add_argument("--fits", "-f", help="fits image to fit")
parser.add_argument(
    "--coords", "-c", help="file containing source pixel " + "coordinates"
)
parser.add_argument("--outfile", "-o", help="file prefix to output sources to")
parser.add_argument(
    "--size",
    "-s",
    type=int,
    help="size of region around " + "each source to include in fit",
    default=25,
)
parser.add_argument(
    "--askap",
    "-a",
    help="file to output ASKAP source " + "positions into for offset calc",
)
args = parser.parse_args()

sources = []
for line in open(args.coords):
    print(line)
    sources.append(tuple(map(int, line[:-1].split(","))))

ex = "jmfitfromfile.py"
sourcefiles = []
size = args.size

for i, s in enumerate(sources):
    print((i, s))
    outfile = args.outfile + f"_source{i}" + ".stats"
    sourcefiles.append(outfile)
    cmd = "{} {} {} {},{},{},{}".format(
        ex,
        args.fits,
        outfile,
        s[0] - size,
        s[1] - size,
        s[0] + size,
        s[1] + size,
    )
    os.system(cmd)


def parse_sourcefile(filename):
    filedict = {}
    for line in open(filename):
        # Lines are formatted like:
        #   FIELD:_VALUE
        # where _ represents at least one space. The values can contain :
        # characters, so we split the lines with ': '
        field, value = tuple(line.split(": "))
        filedict[field] = value[:-1].strip()  # remove newline and whitespace

    return filedict


sourcedicts = []
for sourcefile in sourcefiles:
    sourcedicts.append(parse_sourcefile(sourcefile))

# Save ASKAP source positions
with open(args.askap, "w") as f:
    basestr = "{0},{1},{2},{3}\n"
    for sourcedict in sourcedicts:
        RA = sourcedict["Actual RA"]
        Dec = sourcedict["Actual Dec"]
        # Errors are in mas, convert to arcsec
        RA_err_as = float(sourcedict["Est. RA error (mas)"]) / 1e3
        Dec_err_as = float(sourcedict["Est. Dec error (mas)"]) / 1e3

        f.write(basestr.format(RA, RA_err_as, Dec, Dec_err_as))
