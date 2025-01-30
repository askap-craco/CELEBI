import astropy.io.fits as fits
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord

allinput = input().split(",")
target = allinput[0]
posfileinput = allinput[1]
posfileoutput = allinput[2]

hdulist = fits.open(target)
w = wcs.WCS(hdulist[0].header, hdulist)
hdulist.close()

pos_pix = "#ra,dec\n"

file1 = open(posfileinput, 'r')
lines = file1.readlines()
  
# Strips the newline character
for line in lines:
    pos = line.strip()
    pos = pos.split(',')
    c = SkyCoord(pos[0], pos[1], frame='fk5', unit='deg')
    rapix, decpix = c.to_pixel(wcs=w)
    if c != None:
        new_pos_pix = str(int(rapix.round())) + "," + str(int(decpix.round())) + "\n"
        pos_pix += new_pos_pix

# write to file
with open(posfileoutput, "w") as f:
    f.write(pos_pix)
