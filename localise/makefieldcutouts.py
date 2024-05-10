#this file is for making cutout images of field sources, for independent verification

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

from pathlib import Path
Path("cutouts").mkdir(parents=True, exist_ok=True)

#open the file
filename = 'field.fits'
hdu = fits.open(filename)[0]
data = hdu.data
hdu.header['TIMESYS'] = 'utc' #necessary! it is 'UTC' by default here and this is disagreeable apparently
wcs1 = WCS(hdu.header)
bmaj = hdu.header['BMAJ']*3600 #arcsec
bmin = hdu.header['BMIN']*3600 #arcsec
bpa = hdu.header['BPA'] #deg
pix = np.abs(hdu.header['CDELT1']*3600) #arcsec

#read in field sources by pixel position
fsources = np.loadtxt('field_sources.txt',skiprows=1,delimiter=',')

#size of image
size = (50,50)

#get S/N values
sntxt = open('sources.reg','r')#np.genfromtxt('sources.reg',dtype='str',delimiter='#')
SNlist = []
for line in sntxt.readlines():
    SNlist.append(line.split('S/N=')[1].split('"')[0])
sntxt.close()

#loop over all 50 sources
print('Starting individual field source imaging...')
for i in range(0,50):

    position = (fsources[i][0],fsources[i][1]) #pixel number
    wcs1 = WCS(hdu.header)
    cutout = Cutout2D(data[0,0],position,size,wcs = wcs1.celestial)

    hduCO = hdu.copy()#.header
    hduCO.header.update(cutout.wcs.to_header())
    wcs2 = WCS(hduCO.header) #so we get the coords right for the cutout

    ax = plt.subplot(projection=wcs2,slices=['x', 'y', 0,0])
    ax.imshow(cutout.data, origin='lower')
    #add_beam(ax, major= hdu.header['BMAJ'] * u.deg, minor=hdu.header['BMIN'] * u.deg,
    #         angle=hdu.header['BPA'],hatch='/')#, frame=True)
    beampos = (5,5)
    beam1 = Ellipse(xy=beampos,height=bmaj/pix, width=bmin/pix,
        angle= bpa, edgecolor='orange',facecolor='None',hatch='/') #pixel! because else it is wrong
    ax.add_patch(beam1)
    ax.set_xlabel('RA',fontsize=20)
    ax.set_ylabel('Declination',fontsize=20)
    ax.set_title(f'FieldSource{i}',fontsize=22)
    ax.text(30,45,'S/N = {:0.2f}'.format(float(SNlist[i])),fontsize=16,color='white')
    plt.savefig(f'cutouts/FieldSource{i}.png',bbox_inches='tight')
    plt.close()
print('Field source cutouts images created in cutouts/!')

