#   This script returns the cross matching radius
#   
#                       AB 19 Feb 2025

import sys
import numpy as np
from astropy.io import fits

infitsname  =   sys.argv[1]
matradpar   =   float(sys.argv[2])

with fits.open(infitsname) as filehandle:
    bmaj = filehandle[0].header['BMAJ']
    bmin = filehandle[0].header['BMIN']

radarcsec = matradpar

if (radarcsec < 1.0):
    radarcsec = matradpar*(bmaj*3600.0)

print("Matching radius in arcsec = ", radarcsec)

np.savetxt("match_radius_arcsec.txt", np.array([radarcsec]), fmt="%f")

