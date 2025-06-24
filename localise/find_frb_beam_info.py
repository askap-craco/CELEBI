#####################################################################################

# This script finds beam information of a given FRB and the angular offset of the 
# event w.r.t. the beam centre

#####################################################################################

import argparse
import os
import sys
import re

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord as sc

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 8

#####################################################################################

# SET UP ARGUEMENTS

#####################################################################################

parser = argparse.ArgumentParser(
    description="Find beam information of a given FRB and the angular offset of the event w.r.t. the beam centre"
)
parser.add_argument(
    "--posfile",
    type=str,
    default=None,
    help="FRB position JMFIT file",
)
parser.add_argument(
    "--refvcraftfile",
    type=str,
    default=None,
    help="Reference vcraft file to get beam info",
)
parser.add_argument(
    "--fieldfits",
    type=str,
    default=None,
    help="The FITS image of the field",
)
parser.add_argument(
    "--cenfreqmhz",
    type=float,
    default=1000.0,
    help="Central frequency in MHz",
)
parser.add_argument(
    "--bwmhz",
    type=float,
    default=336.0,
    help="Band width in MHz",
)
parser.add_argument(
    "--diametre",
    type=float,
    default=12.0,
    help="Dish diameter in metre",
)
parser.add_argument(
    "--plotpng",
    type=str,
    default=None,
    help="Name of the plot file without extension",
)

args = parser.parse_args()
#print(f"\nArguments specified: {args}\n")

#####################################################################################

def get_adj_beams(bid, footprint):
    
    # Probably the dumbest function AB has ever written   
    
    adjbeams    =   0
    position    =   "Unknown footprint! Contact Ryan M Shannon."

    if (footprint == "closepack36"):
        if (bid in [0,35]):
            adjbeams = 2
            position = "corner"
        elif (bid in [5,11,12,23,24,30]):
            adjbeams = 3
            position = "semi-corner"
        elif (bid in [1,2,3,4,31,32,33,34]):
            adjbeams = 4
            position = "edge"
        elif (bid in [6,17,18,29]):
            adjbeams = 5
            position = "near-edge"
        else:
            adjbeams = 6
            position = "good"

    elif (footprint == "square_6x6"):
        if (bid<=16):
            adjbeams = 4
            position = "good"
        elif (bid in [16, 21, 26, 31]):
            adjbeams = 2
            position = "corner"
        else:
            adjbeams = 3
            position = "edge"

    return (adjbeams,position)

#   --------------------------------------------------------------------------------

def stringToRad(posstr, was_hms):
    posstr = posstr.strip()
    if was_hms and "h" in posstr:
        posstr = posstr.replace("h", ":")
        posstr = posstr.replace("m", ":")
        posstr = posstr.replace("s", "")
    if not was_hms and "d" in posstr:
        posstr = posstr.replace("d", ":")
        posstr = posstr.replace("'", ":")
        posstr = posstr.replace("\"", "")
    splits = posstr.split(':')
    try:
        d = int(splits[0])
        m = 0
        s = 0
        negmult = 1.0
        if len(splits) > 1 and len(splits[1].strip()) > 0:
            m = int(splits[1])
        if len(splits) > 2 and len(splits[2].strip()) > 0:
            s = float(splits[2])
        if posstr[0] == "-":
            d = -d
            negmult = -1.0
        radval = d + m/60.0 + s/3600.0
        radval *= np.pi/180.0
        radval *= negmult
        if was_hms:
            radval *= 180.0/12.0
    except ValueError:
        print("Bad position string", posstr)
        radval = -999
    return(radval)

#   -------------------------------------------------------------------------------

def get_sources_pos(srcposfile):

    fields = {}
    f = open(srcposfile)

    for line in f:
        line = line[:-1].split(": ")  # [:-1] to trim newline
        fields[line[0]] = line[1].strip()  # strip extra whitespace
    f.close()

    ra_hms = fields["Actual RA"]  # hms
    dec_dms = fields["Actual Dec"]  # dms

    rarad = stringToRad(ra_hms, True)
    decrad = stringToRad(dec_dms, False)

    return (rarad, decrad)

#   --------------------------------------------------------------------------------

def posdiff(targetrarad, targetdecrad, calrarad, caldecrad):
    sinsqdecdiff = np.sin((targetdecrad-caldecrad)/2.0)
    sinsqdecdiff = sinsqdecdiff*sinsqdecdiff
    sinsqradiff  = np.sin((targetrarad-calrarad)/2.0)
    sinsqradiff  = sinsqradiff*sinsqradiff

    return(2*np.arcsin(np.sqrt(sinsqdecdiff + np.cos(targetdecrad)*np.cos(caldecrad)*sinsqradiff)))

#   --------------------------------------------------------------------------------

def get_field_centre(fieldfitsfile):
    fldfits  = fits.open(fieldfitsfile)
    fldhdr   = fldfits[0].header
    fldra    = np.deg2rad( fldhdr['CRVAL1'] - (fldhdr['CRPIX1'] - fldhdr['NAXIS1']*0.5) * fldhdr['CDELT1'] )   
    flddec   = np.deg2rad( fldhdr['CRVAL2'] - (fldhdr['CRPIX2'] - fldhdr['NAXIS2']*0.5) * fldhdr['CDELT2'] ) 
    fldfits.close()
    return (fldra, flddec)

#   --------------------------------------------------------------------------------

class VcraftHdr:
    def __init__(self, hdrfpath):
        self.hdrfpath = hdrfpath
        self.__init_load()
        self._load_hdrcontent()

    def __init_load(self):
        self.hdrcontent = self._load_hdr_content()
        self.frbname = self._load_frb_name()


    def _load_hdr_content(self, ):
        with open(self.hdrfpath, ) as fp:
            content = fp.read()
        return content

    def _load_frb_name(self,):
        return self.hdrfpath.split("/")[4]

    ### load information from header
    def _load_beam_info(self):
        match = re.findall(r"BEAM (\d+) #", self.hdrcontent)
        if len(match) == 0: return -1
        return int(match[0])

    def _load_footprint(self):
        match = re.findall("FOOTPRINT (.+) #", self.hdrcontent)
        if len(match) == 0: return ",,"
        return match[0]

    def _load_beam_pointing(self):
        ### ra
        match = re.findall("BEAM_RA (.+) #", self.hdrcontent)
        ra = 999. if len(match) == 0 else float(match[0])
        ### dec
        match = re.findall("BEAM_DEC (.+) #", self.hdrcontent)
        dec = 999. if len(match) == 0 else float(match[0])
        return ra, dec

    def _load_hdrcontent(self):
        self.vcraftbeam = self._load_beam_info()
        ### get footprint
        self.footprint, self.pitch, self.pafang = self._load_footprint().split(",")
        self.beamra, self.beamdec = self._load_beam_pointing()
    
#   --------------------------------------------------------------------------------

fieldradec  =   get_field_centre(args.fieldfits)
frbradec    =   get_sources_pos(args.posfile)
offindeg    =   np.rad2deg(posdiff(frbradec[0],frbradec[1],fieldradec[0],fieldradec[1]))

vcrafthdr   =   VcraftHdr(args.refvcraftfile)
adjbms,bpos =   get_adj_beams(int(vcrafthdr.vcraftbeam) // 2, vcrafthdr.footprint)

fmhzarr     =   np.arange(args.cenfreqmhz - (args.bwmhz/2), args.cenfreqmhz + (args.bwmhz/2) + 0.5, 1.0, dtype=float)
nulldeg     =   np.rad2deg(1.22 * (299.792458/fmhzarr) / args.diametre)
hppdeg      =   np.rad2deg(0.5 * (299.792458/fmhzarr) / args.diametre)

fig 	= plt.figure(figsize=(4.0,3.6))
ax 		= fig.add_axes([0.12, 0.10, 0.85,0.82])
ax.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)

ax.set_title("Beam "+str(int(vcrafthdr.vcraftbeam) // 2)+" in "+vcrafthdr.footprint+" [ "+bpos+" beam ]")
ax.plot(fmhzarr/1.0e3, nulldeg, 'r--', label='First null')
ax.plot(fmhzarr/1.0e3, hppdeg, 'c-', label='Half power point')
ax.axhline(y=offindeg,c='k',ls=':',label="FRB position")

ax.set_ylim(ymin=0)
ax.set_xlim([args.cenfreqmhz*1.0e-3 - (args.bwmhz/2)*1.0e-3, args.cenfreqmhz*1.0e-3 + (args.bwmhz/2)*1.0e-3])
ax.legend(loc='lower left')
ax.set_ylabel(r'Angular distance from centre (deg)')
ax.set_xlabel(r'Frequency (GHz)')
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_label_coords(-0.09, 0.5)

plt.savefig(args.plotpng+".png", transparent=True, format='png')
plt.close()

print(f"\n    ASKAP Beam information\n")
print(f"Beam ID       ", int(vcrafthdr.vcraftbeam) // 2)
print(f"Footprint     ", vcrafthdr.footprint)
print(f"Pitch         ", vcrafthdr.pitch)
print(f"PAF angle     ", vcrafthdr.pafang)
print(f"Beam Centre   ", vcrafthdr.beamra, vcrafthdr.beamdec, "(for an arbitrary antenna) \n")
print(f"Adjacent beams", adjbms, "["+bpos+"]")

print(f"\n    Relative FRB position\n")
print(f"Field centre  ", np.rad2deg(fieldradec), "deg")
print(f"FRB offset    ", offindeg," deg from the field centre")
print(f"See figure    ", args.plotpng, ".png for relative position")





