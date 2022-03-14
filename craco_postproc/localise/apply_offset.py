#!/usr/bin/env python3

###############################################################################
# craftpy2 stage 2.6: apply offset
# Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
# Last modified: 2021-05-24
#
# Applies the systematic offset calculated from field source positions and
# applies it to the FRB position
###############################################################################

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from astropy.coordinates import SkyCoord as sc
from astropy import units as un
import numpy as np

parser = ArgumentParser(description='Apply the systematic offset calculated '
                                    'from field source positions and apply it '
                                    'to the FRB position.',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('--frb', help='Pre-offset FRB JMFIT stats file')
parser.add_argument('--offset', help='Offset file')
args = parser.parse_args()

class Coord(object):
    def __init__(self, stats):
        fields = {}
        f = open(stats)
        for line in f:
            line = line[:-1].split(': ')        # [:-1] to trim newline
            fields[line[0]] = line[1].strip()   # strip extra whitespace
        f.close()

        self.ra_hms = fields['Actual RA']                       # hms
        self.ra_err = float(fields['Est. RA error (mas)'])/1e3  # arcseconds
        self.dec_dms = fields['Actual Dec']                     # dms
        self.dec_err = float(fields['Est. Dec error (mas)'])/1e3# arcseconds

FRB = Coord(args.frb)
offset_ra, offset_ra_err, offset_dec, offset_dec_err = np.loadtxt(args.offset)

frb_pos = sc(FRB.ra_hms, FRB.dec_dms, unit=('hourangle', 'deg'))
frb_err = sc(FRB.ra_err, FRB.dec_err, unit='arcsecond')

offset = sc(offset_ra, offset_dec, unit='arcsecond')
offset_err = sc(offset_ra_err, offset_dec_err, unit='arcsecond')

# sum positions
final_pos = sc(frb_pos.ra + offset.ra, frb_pos.dec, offset.dec)
final_err = sc((frb_err.ra**2 + offset_err.ra**2)**0.5,
               (frb_err.dec**2 + offset_err.dec**2)**0.5)

print('FINAL FRB POSITION')
print(f'RA:\t{final_pos.ra.to_string("hour")}'
      f'\t+-{final_err.ra.to_string("arcsecond")}')
print(f'Dec:\t{final_pos.dec.to_string("deg")}'
      f'\t+-{final_err.dec.to_string("arcsecond")}')
