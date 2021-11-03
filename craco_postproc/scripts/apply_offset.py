#!/usr/bin/env python3

###############################################################################
# craftpy2 stage 2.6: apply offset
# Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
# Last modified: 2021-05-24
#
# Applies the systematic offset calculated from field source positions and
# applies it to the FRB position
###############################################################################

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from astropy import units as un
from astropy.coordinates import SkyCoord as sc

parser = ArgumentParser(
    description="Apply the systematic offset calculated "
    "from field source positions and apply it "
    "to the FRB position.",
    formatter_class=ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--ra_frb", help="Pre-offset FRB right ascension (hms)")
parser.add_argument("--dec_frb", help="Pre-offset FRB declination (dms)")
parser.add_argument(
    "--ra_err_frb", help="Uncertainty in FRB right ascension " "(arcsec)"
)
parser.add_argument(
    "--dec_err_frb", help="Uncertainty in FRB declination " "(arcsec)"
)
parser.add_argument("--ra_offset", help="Offset in right ascension (arcsec)")
parser.add_argument("--dec_offset", help="Offset in declination (arcsec)")
parser.add_argument(
    "--ra_err_offset",
    help="Uncertainty in offset in right " "ascension (arcsec)",
)
parser.add_argument(
    "--dec_err_offset", help="Uncertainty in offset in " "declination (arcsec)"
)
args = parser.parse_args()

frb_pos = sc(args.ra_frb, args.dec_frb, unit=("hourangle", "deg"))
frb_err = sc(args.ra_err_frb, args.dec_err_frb, unit="arcsecond")

offset = sc(args.ra_offset, args.dec_offset, unit="arcsecond")
offset_err = sc(args.ra_err_offset, args.dec_err_offset, unit="arcsecond")

# sum positions
final_pos = sc(frb_pos.ra + offset.ra, frb_pos.dec, offset.dec)
final_err = sc(
    (frb_err.ra ** 2 + offset_err.ra ** 2) ** 0.5,
    (frb_err.dec ** 2 + offset_err.dec ** 2) ** 0.5,
)

print("FINAL FRB POSITION")
print(
    f'RA:\t{final_pos.ra.to_string("hour")}'
    f'\t+-{final_err.ra.to_string("arcsecond")}'
)
print(
    f'Dec:\t{final_pos.dec.to_string("deg")}'
    f'\t+-{final_err.dec.to_string("arcsecond")}'
)
