import sys
from astropy.coordinates import SkyCoord as sc

ra_deg = sys.argv[1]
dec_deg = sys.argv[2]

radec = sc(ra_deg, dec_deg, unit='deg')

ra_hms = radec.ra.to_string('hour')
dec_dms = radec.dec.to_string('deg')

print('%s %s' % (ra_hms, dec_dms))