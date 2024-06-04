from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

parser = ArgumentParser(
    description="Figure out the galactic coordinates of the FRB and choose systematic uncertainties",
    formatter_class=ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--frbra", help="FRB right ascension string")
parser.add_argument("--frbdec", help="FRB declination string")
parser.add_argument("--planelatcut", type=float, help="Galactic latitude defining 'on-plane' region")
parser.add_argument("--onplanera", help="RA systematic on-plane")
parser.add_argument("--onplanedec", help="RA systematic off-plane")
parser.add_argument("--offplanera", help="Dec systematic on-plane")
parser.add_argument("--offplanedec", help="Dec systematic off-plane")

args = parser.parse_args()

frb = SkyCoord(args.frbra, args.frbdec, unit=(u.hourangle, u.deg), frame='icrs')
if np.fabs((frb.galactic.b * u.deg).value) < args.planelatcut:
    print(args.onplanera, args.onplanedec)
else:
    print(args.offplanera, args.offplanedec)
