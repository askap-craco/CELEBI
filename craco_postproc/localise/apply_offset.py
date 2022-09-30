###############################################################################
# craftpy2 stage 2.6: apply offset
# Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
# Last modified: 2021-05-24
#
# Applies the systematic offset calculated from field source positions and
# applies it to the FRB position
###############################################################################

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np
from astropy import units as un
from astropy.coordinates import SkyCoord as sc

parser = ArgumentParser(
    description="Apply the systematic offset calculated "
    "from field source positions and apply it "
    "to the FRB position.",
    formatter_class=ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--frb", help="Pre-offset FRB JMFIT stats file")
parser.add_argument("--offset", help="Offset file")
args = parser.parse_args()


class Coord:
    def __init__(self, stats):
        fields = {}
        f = open(stats)
        for line in f:
            line = line[:-1].split(": ")  # [:-1] to trim newline
            fields[line[0]] = line[1].strip()  # strip extra whitespace
        f.close()

        self.ra_hms = fields["Actual RA"]  # hms
        self.ra_err = float(fields["Est. RA error (mas)"]) / 1e3  # arcseconds
        self.dec_dms = fields["Actual Dec"]  # dms
        self.dec_err = (
            float(fields["Est. Dec error (mas)"]) / 1e3
        )  # arcseconds


FRB = Coord(args.frb)
offset_ra, offset_ra_err, offset_dec, offset_dec_err = np.loadtxt(args.offset)

frb_pos = sc(FRB.ra_hms, FRB.dec_dms, unit=("hourangle", "deg"))

final_pos = sc(
    frb_pos.ra + offset_ra*un.arcsecond/np.cos(frb_pos.dec.to(un.rad)),
    frb_pos.dec + offset_dec*un.arcsecond
)

final_err_ra = (FRB.ra_err ** 2 + offset_ra_err ** 2) ** 0.5
final_err_dec = (FRB.dec_err ** 2 + offset_dec_err ** 2) ** 0.5


print("FINAL FRB POSITION")
print(
    f'RA:\t{final_pos.ra.to_string("hour")}'
    f'\t+-{final_err_ra}'
)
print(
    f'Dec:\t{final_pos.dec.to_string("deg")}'
    f'\t+-{final_err_dec}'
)
