#this file is for making cutout images of field sources, for independent verification

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import sys

from pathlib import Path

# Read inputs
filename = sys.argv[1]
dobm     = int(sys.argv[2])

#open the file
hdu = fits.open(filename+'.fits')[0]
data = hdu.data
hdu.header['TIMESYS'] = 'utc' #necessary! it is 'UTC' by default here and this is disagreeable apparently
wcs1 = WCS(hdu.header)
pix = np.abs(hdu.header['CDELT1']*3600) #arcsec

fig = plt.figure(figsize=(12,12),dpi=300)
ax = plt.subplot(projection=wcs1,slices=['x', 'y', 0,0])
ax.imshow(data[0,0], origin='lower')

if (dobm):
    bmaj = hdu.header['BMAJ']*3600 #arcsec
    bmin = hdu.header['BMIN']*3600 #arcsec
    bpa = hdu.header['BPA'] #deg
    beampos = (20,20)
    beam1 = Ellipse(xy=beampos,height=bmaj/pix, width=bmin/pix,
        angle= bpa, edgecolor='w',facecolor='None') #pixel! because else it is wrong
    ax.add_patch(beam1)
ax.set_xlabel('RA',fontsize=20)
ax.set_ylabel('Declination',fontsize=20)
plt.savefig(filename+'.png',bbox_inches='tight')
plt.close()


