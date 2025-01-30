#   This script ALMAfies the ASKAP visibilities
#   by just renaming ASKAP as ALMA !!
#   
#   Written by AB based on codes by MG, Dec 1, 2023

import sys
from astropy.io import fits

infitsname  =   sys.argv[1]

with fits.open(infitsname, mode='update') as filehandle:
    filehandle[0].header['INSTRUME'] = 'ALMA'
    filehandle[0].header['TELESCOP'] = 'ALMA'
