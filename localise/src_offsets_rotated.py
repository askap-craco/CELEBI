#####################################################################################

# GENERAL PYTHON IMPORTS

#####################################################################################

import argparse
# import glob
import os
# import re
import sys
import math
from itertools import dropwhile #, islice

import matplotlib.pyplot as plt
import numpy as np
# from astropy import constants as const
# from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord as sc

#####################################################################################

# SET UP ARGUEMENTS

#####################################################################################

# Set up parser with description string
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""\
PURPOSE: Continuum source position offset calculator and plotter

This script calculates the offsets between continuum field source positions from ASKAP and 
one or more source catalogues and their uncertainties.

INPUT:
File: containing the source positions from ASKAP
File: containing the ASKAP source names
File: containing relevant jmfit filenames
File(s): containing catalogue source positions for comparison
Notes on file formats:
txt files with comma separated values as follows
ASKAP: RA(hms),RA_err(ms),Dec(dms),Dec_err(mas)
NVSS: RA(hms),RA_err(ms),Dec(dms),Dec_err(arcsec)
SUMSS: RA(hms),RA_err(arcsec),Dec(dms),Dec_err(arcsec)
FIRST: RA(hms),Dec(dms),Fpeak(mJy/bm),rms(mJy/bm),Maj("),Min("),PA(deg),fMaj("),fMin(")
VLASS: RA(hms),RA_err(ms),Dec(dms),Dec_err(mas)

OUTPUT:
Plot: source offsets (ASKAP vs. Catalogue) with weighted mean
File: .dat file containing the source offsets and uncertainties to be used for
weighted_multi_image_fit.py infile1 infile2 ...
  e.g. python multi_image_fit.py 190102.dat 190608.dat 190711.dat
      Infiles have: one column per source, rows are ra_offset, dec_offset, ra_err, dec_err
      Please use simple text file with numbers and whitespace separation

Revision history:
     Written by Cherie Day on 5 June, 2020
     Python version: 3.7.0

     Modified by AB on March 14, 2023
""",
)
parser.add_argument(
    "--askappos",
    type=str,
    default=None,
    help="The file containing the ASKAP continuum source positions",
)
parser.add_argument(
    "--askapnames",
    type=str,
    default=None,
    help="The file containing the ASKAP continuum source names",
)
parser.add_argument(
    "--jmfitnames",
    type=str,
    default=None,
    help="The file containing the list of jmfit files",
)
parser.add_argument(
    "--fieldfits",
    type=str,
    default=None,
    help="The FITS image of the field",
)
parser.add_argument(
    "--first",
    type=str,
    default=None,
    help="The file containing the FIRST catalogue sources if use",
)
parser.add_argument(
    "--nvss",
    type=str,
    default=None,
    help="Set if the NVSS catalogue is being used",
)
parser.add_argument(
    "--racs",
    type=str,
    default=None,
    help="Set if the RACS catalogue is being used",
)
parser.add_argument(
    "--sumss",
    type=str,
    default=None,
    help="Set if the SUMSS catalogue is being used",
)
parser.add_argument(
    "--vlass",
    type=str,
    default=None,
    help="Set if the VLASS catalogue is being used",
)
parser.add_argument(
    "--frbtitletext",
    type=str,
    default=None,
    help="The name of the ASKAP FRB field being compared",
)

args = parser.parse_args()
print(f"\nArguments specified: {args}\n")

#####################################################################################

# CATCH ERRORS

#####################################################################################

if len(sys.argv) < 2:
    parser.print_usage()
    sys.exit()

if args.askappos is None:
    parser.error(
        "You must specify a text file containing the ASKAP source information"
    )

if args.askapnames is None:
    parser.error("You must specify a text file containing the ASKAP names")

if (
    (args.first is None)
    and (args.nvss is None)
    and (args.sumss is None)
    and (args.vlass is None)
    and (args.racs is None)
):
    parser.error(
        "You must specify one or more text files containing the catalogue source information"
    )

if args.frbtitletext is None:
    parser.error("You must specify a title for the plots")


#####################################################################################

# GLOBAL PARAMETER DEFINITIONS

#####################################################################################

# Get file name variables for ASKAP file and any catalogue files selected
askap_posfile = str(args.askappos)
if args.first is not None:
    first_posfile = str(args.first)
if args.nvss is not None:
    nvss_posfile = str(args.nvss)
if args.racs is not None:
    racs_posfile = str(args.racs)
if args.sumss is not None:
    sumss_posfile = str(args.sumss)
if args.vlass is not None:
    vlass_posfile = str(args.vlass)

np.set_printoptions(precision=15)

# Set up colour dictionary with colourblind friendly scheme
col_bank = {
    "1": "#40004b",
    "2": "#543005",
    "3": "#9970ab",
    "4": "#bf812d",
    "5": "#8c510a",
    "6": "#003c30",
    "7": "#a6dba0",
    "8": "#5aae61",
    "9": "#1b7837",
    "10": "#00441b",
    "11": "#af8dc3",
    "12": "#7fbf7b",
    "13": "#008837",
    "14": "#7b3294",
    "15": "#a6611a",
    "16": "#018571",
}

RAD_IN_ARCSEC = 180.0*3600/np.pi

#####################################################################################

# DEFINE FUNCTIONS USED FOR OFFSET AND UNCERTAINTY CALCULATIONS

#####################################################################################


def is_comment(s):
    """
    function to check is a line is a comment

    Parameters
    ----------
    s : str
        Line in a file being read

    Returns
    ----------
    True
        if a line starts with #
    """

    return s.startswith("#")


def grabsrcinfo(infofile):

    """
    This function reads a file containing source information
    and returns each line as a list entry

    Parameters
    ----------
    infofile : str
        File containing source information

    Returns
    ----------
    list
        a list of strings
    """

    source_info = []

    with open(infofile) as info:
        # Ignore comments and add each line into a list
        for currentline in dropwhile(is_comment, info):
            source_info.append(currentline)

    return source_info


def getradec(src_info, first=False):

    """
    This function takes a list of source information,
    grabs the positions, and converts them to arcseconds.

    Parameters
    ----------
    src_info : list
        List containing source information
    first : boolean
        True if FIRST catalogue; default False

    Returns
    ----------
    dict of lists
        ra_arcsec : list with the RAs for each source in arcsecond
        dec_arcsec : list with the Decs for each source in arcsecond
        dec_radian : list of Decs for each source in radians

    """
    pos_dict = {}

    # Extract positions in hms and dms
    ra_hms = [src_info[i].split(",")[0] for i in np.arange(len(src_info))]
    print(ra_hms)
    try:
        ra_hms = [
            (
                ra.split(":")[0]
                + "h"
                + ra.split(":")[1]
                + "m"
                + ra.split(":")[2]
                + "s"
            )
            for ra in ra_hms
        ]
    except:
        ra_hms = ra_hms
    print(ra_hms)

    if first is True:
        dec_dms = [src_info[i].split(",")[1] for i in np.arange(len(src_info))]
    else:
        dec_dms = [src_info[i].split(",")[2] for i in np.arange(len(src_info))]

    try:
        dec_dms = [
            (
                dec.split(":")[0]
                + "d"
                + dec.split(":")[1]
                + "m"
                + dec.split(":")[2]
                + "s"
            )
            for dec in dec_dms
        ]
    except:
        dec_dms = dec_dms

    # Convert RA, Dec into arcsec
    ra_arcsec = [
        (sc(ra, dec, frame="icrs").ra.arcsec)
        for ra, dec in zip(ra_hms, dec_dms)
    ]
    dec_arcsec = [
        (sc(ra, dec, frame="icrs").dec.arcsec)
        for ra, dec in zip(ra_hms, dec_dms)
    ]

    # Convert Dec to radians
    dec_rad = [
        (sc(ra, dec, frame="icrs").dec.radian)
        for ra, dec in zip(ra_hms, dec_dms)
    ]

    # Add positions to dictionary
    pos_dict["ra_arc"] = ra_arcsec
    pos_dict["dec_arc"] = dec_arcsec
    pos_dict["dec_rad"] = dec_rad

    return pos_dict


def getfirst_unc(src_info):

    """
    This function takes a list of FIRST source information
    and determines the 90% uncertainty in RA and Dec using the
    following empirical expressions from the FIRST catalogue:

    unc(90% confidence) = Size  (1/SNR + 1/20)  arcsec
    SNR = (Fpeak-0.25) / RMS

    It then converts these uncertainties to 60% to match
    the other catalogue and fitted uncertainties

    Parameters
    ----------
    src_info : list
        List containing source information

    Returns
    ----------
    dict of lists
        ra_unc_arcsec : list with the RA uncertainty for each source in arcsecond
        dec_unc_arcsec : list with the Dec uncertainty for each source in arcsecond

    """

    first_unc_dict = {}

    # Grab all the relevant values from the source info
    # Peak flux density: mJy/beam
    pflux = [
        #np.
        float(src_info[i].split(",")[2]) for i in np.arange(len(src_info))
    ]
    # Flux denisity RMS: mJy/beam
    rms = [
        #np.
        float(src_info[i].split(",")[3]) for i in np.arange(len(src_info))
    ]
    # Measured major axis: arcsec
    fmaj = [
        #np.
        float(src_info[i].split(",")[7]) for i in np.arange(len(src_info))
    ]
    # Measured minor axis: arcsec
    fmin = [
        #np.
        float(src_info[i].split(",")[8]) for i in np.arange(len(src_info))
    ]
    # Measured position angle: degrees
    fpa = [
        #np.
        float(src_info[i].split(",")[9]) for i in np.arange(len(src_info))
    ]

    # Project FWHM major and minor axes onto RA, Dec axes
    ra_proj = [
        np.sqrt(
            (np.cos(fpa[i] * np.pi / 180) ** 2 / fmaj[i] ** 2)
            + (np.sin(fpa[i] * np.pi / 180) ** 2 / fmin[i] ** 2)
        )
        ** (-1)
        for i in np.arange(len(fpa))
    ]
    dec_proj = [
        np.sqrt(
            (np.sin(fpa[i] * np.pi / 180) ** 2 / fmaj[i] ** 2)
            + (np.cos(fpa[i] * np.pi / 180) ** 2 / fmin[i] ** 2)
        )
        ** (-1)
        for i in np.arange(len(fpa))
    ]

    # Calculate the SNR
    snr = [(pflux[i] - 0.25) / rms[i] for i in np.arange(len(pflux))]

    # Derive the 90% confidence uncertainties for RA and Dec in arcsec
    ra_unc_arcsec_90 = [
        ra_proj[i] * ((1 / snr[i]) + (1 / 20.0))
        for i in np.arange(len(ra_proj))
    ]
    dec_unc_arcsec_90 = [
        dec_proj[i] * ((1 / snr[i]) + (1 / 20.0))
        for i in np.arange(len(dec_proj))
    ]

    # Convert to 60% confidence uncertainties for RA and Dec in arcsec
    ra_unc_arcsec = np.array(ra_unc_arcsec_90) / 1.645
    dec_unc_arcsec = np.array(dec_unc_arcsec_90) / 1.645

    first_unc_dict["ra_unc_arcsec"] = ra_unc_arcsec
    first_unc_dict["dec_unc_arcsec"] = dec_unc_arcsec

    return first_unc_dict


def getradec_unc(
    src_info,
    dec_rad,
    askap=False,
    first=False,
    nvss=False,
    sumss=False,
    vlass=False,
    racs=False,
):

    """
    This function takes a list of source information,
    grabs the positional uncertainties, and converts them to arcseconds.

    If first=True, calls getfirst_unc(src_info) to derive the RA, Dec 90%
    uncertainties.

    Parameters
    ----------
    src_info : list
        List containing source information
    dec_rad : list
        List of Declinations (floats) in radians
    askap : boolean
        True if ASKAP source; default False
    first : boolean
        True if FIRST catalogue; default False
    nvss : boolean
        True if NVSS catalogue; default False
    sumss : boolean
        True if SUMSS catalogue; default False
    vlass : boolean
        True if VLASS sources; default False

    Returns
    ----------
    dict of lists
        ra_unc_arcsec : list with the RA uncertainties for each source in arcsecond
        dec_unc_arcsec : list with the Dec uncertainties for each source in arcsecond

    """

    unc_dict = {}

    # Extract the uncertainties in the default unit of the source info and convert to arcsec if needed
    # ASKAP or NVSS
    if (askap is True) or (nvss is True) or (vlass is True) or (racs is True):
        # TEMP: RA_err expected to be arcseconds
        ra_unc_arcsec = [
            #np.
            float(src_info[i].split(",")[1])
            for i in np.arange(len(src_info))
        ]

        """
        # NVSS default RA_err unit: s
        ra_unc_s = [np.float(src_info[i].split(',')[1]) for i in np.arange(len(src_info))]
        print(ra_unc_s)
        # Convert to arcsec: s * 15cos(dec_rad)
        ra_unc_arcsec = [ra_unc_s[i] * 15 * np.cos(dec_rad[i]) for i in np.arange(len(ra_unc_s))]
        print(15*np.cos(dec_rad))
        print(ra_unc_arcsec)
        """

        ### Old: RA_err default unit: ms
        ### ra_unc_ms = [np.float(src_info[i].split(',')[1]) for i in np.arange(len(src_info))]
        ### Convert to arcsec: ms * 15cos(dec_rad) / 1000
        ### ra_unc_arcsec = [ra_unc_ms[i] * 15 * np.cos(dec_rad[i]) / 1000. for i in np.arange(len(ra_unc_ms))]
    if (askap is True) or (vlass is True):
        # Dec_err default unit: mas
        # dec_unc_mas = [np.float(src_info[i].split(',')[3]) for i in np.arange(len(src_info))]
        # dec_unc_arcsec = [dec_unc_mas[i]/1000. for i in np.arange(len(src_info))]
        dec_unc_arcsec = [
            #np.
            float(src_info[i].split(",")[3])
            for i in np.arange(len(src_info))
        ]
    if (nvss is True) or (sumss is True) or (racs is True):
        # Dec_err default unit: arcsec
        dec_unc_arcsec = [
            #np.
            float(src_info[i].split(",")[3])
            for i in np.arange(len(src_info))
        ]
    if first is True:
        ra_unc_arcsec = getfirst_unc(src_info)["ra_unc_arcsec"]
        dec_unc_arcsec = getfirst_unc(src_info)["dec_unc_arcsec"]

    print(ra_unc_arcsec)
    print(dec_unc_arcsec)

    # Add uncertainties to dictionary
    unc_dict["ra_unc_arcsec"] = ra_unc_arcsec
    unc_dict["dec_unc_arcsec"] = dec_unc_arcsec

    return unc_dict


def getoffsets(ra_arcsec1, dec_arcsec1, ra_arcsec2, dec_arcsec2, dec2_rad):

    """
    This function determines the offset between two source positions

    Parameters
    ----------
    ra_arcsec1 : list
        List containing source 1's RA
    dec_arcsec1 : list
        List containing source 1's Dec
    ra_arcsec2 : list
        List containing source 2's RA
    dec_arcsec2 : list
        List containing source 2's Dec
    dec_ra : list
        List containing source 2's Dec in radians

    Returns
    ----------
    dict of lists
        ra_offset_arcsec : list with the RA uncertainties for each source in arcsecond
        dec_offset_arcsec : list with the Dec uncertainties for each source in arcsecond

    """

    offsets_dict = {}

    # Calculate the RA, Dec offsets
    ra_offset_arcsec = [
        (ra2 - ra1) * np.cos(dec)
        for ra2, ra1, dec in zip(ra_arcsec2, ra_arcsec1, dec2_rad)
    ]
    dec_offset_arcsec = [
        dec2 - dec1 for dec2, dec1 in zip(dec_arcsec2, dec_arcsec1)
    ]

    offsets_dict["ra_offset_arcsec"] = ra_offset_arcsec
    offsets_dict["dec_offset_arcsec"] = dec_offset_arcsec

    return offsets_dict


def gettotaloffsetunc(
    ra_unc_arcsec1, dec_unc_arcsec1, ra_unc_arcsec2, dec_unc_arcsec2
):

    """
    This function determines the total uncertainty in the offset between two source positions

    Parameters
    ----------
    ra_unc_arcsec1 : list
        List containing source 1's RA uncertainty in arcsec
    dec_unc_arcsec1 : list
        List containing source 1's Dec uncertainty in arcsec
    ra_unc_arcsec2 : list
        List containing source 2's RA uncertainty in arcsec
    dec_unc_arcsec2 : list
        List containing source 2's Dec uncertainty in arcsec

    Returns
    ----------
    dict of lists
        ra_offset_unc_arcsec : list with the RA uncertainties for each source in arcsecond
        dec_offset_unc_arcsec : list with the Dec uncertainties for each source in arcsecond

    """

    offset_unc_dict = {}

    # Calculate total offset uncertainty
    totra_unc = [
        np.sqrt(ra_unc2 ** 2 + ra_unc1 ** 2)
        for ra_unc2, ra_unc1 in zip(ra_unc_arcsec2, ra_unc_arcsec1)
    ]
    totdec_unc = [
        np.sqrt(dec_unc2 ** 2 + dec_unc1 ** 2)
        for dec_unc2, dec_unc1 in zip(dec_unc_arcsec2, dec_unc_arcsec1)
    ]

    offset_unc_dict["totra_unc"] = totra_unc
    offset_unc_dict["totdec_unc"] = totdec_unc

    return offset_unc_dict


#--------------------------------------------------------------------

# Astro utility functions

#--------------------------------------------------------------------

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
        radval *= math.pi/180.0
        radval *= negmult
        if was_hms:
            radval *= 180.0/12.0
    except ValueError:
        print("Bad position string", posstr)
        radval = -999
    return(radval)


def posdiff(targetrarad, targetdecrad, calrarad, caldecrad):
    sinsqdecdiff = math.sin((targetdecrad-caldecrad)/2.0)
    sinsqdecdiff = sinsqdecdiff*sinsqdecdiff
    sinsqradiff  = math.sin((targetrarad-calrarad)/2.0)
    sinsqradiff  = sinsqradiff*sinsqradiff

    return(2*math.asin(math.sqrt(sinsqdecdiff + math.cos(targetdecrad)*math.cos(caldecrad)*sinsqradiff)))


def posradians2string(rarad, decrad):
    rah = rarad * 12 / math.pi
    rhh = int(rah)
    rmm = int(60*(rah - rhh))
    rss = 3600*rah - (3600*rhh + 60*rmm)
    decd = decrad * 180 / math.pi
    decformat = "+%02d:%02d:%010.7f"
    if decd < 0:
        decd = -decd
        decformat = '-' + decformat[1:]
    ddd = int(decd)
    dmm = int(60*(decd - ddd))
    dss = 3600*decd - (3600*ddd + 60*dmm)
    rastring  = "%02d:%02d:%011.8f" % (rhh,rmm,rss)
    decstring = decformat % (ddd, dmm, dss)
    return(rastring, decstring)



class AstroObsResult:
    def __init__(self, obj, exp, lines, isexp):
        self.object = obj
        self.experiment = exp
        self.isexploratory = isexp
        self.gategain = -1.0
        self.lowerlimit = True
        self.mjd = float(lines[0].split()[3][:-1])
        self.raerrmas = float(lines[10].split()[-1])
        self.raerrhhh = float(lines[11].split()[-1])
        self.decerrmas = float(lines[12].split()[-1])
        self.snr = float(lines[4].split()[-1])
        self.flux = float(lines[3].split()[-1])
        self.ra = lines[7].split()[-1]
        self.dec = lines[8].split()[-1]
        self.rarad = 0
        self.decrad = 0
        if self.snr != 0.0:
            self.rarad = stringToRad(self.ra, True)
            self.decrad = stringToRad(self.dec, False)
        tempsplit = lines[9].split()
        self.fitmin = float(tempsplit[1].split('x')[0])
        self.fitmaj = float(tempsplit[1].split('x')[1])
        self.fitpa  = float(tempsplit[3])
        self.hasbeam = False
        if len(tempsplit) > 6 and "x" in tempsplit[6]:
            self.hasbeam = True
            self.beammin = float(tempsplit[6].split('x')[0])
            self.beammaj = float(tempsplit[6].split('x')[1])
            self.beampa  = float(tempsplit[8])
            self.beamparad = math.pi*self.beampa/180.0

    def updateRA(self, newrarad, newraraderr=-1):
        self.rarad = newrarad
        newra, xxx = posradians2string(newrarad, 1)
        self.ra = newra
        if newraraderr > 0:
            self.raerrmas = newraraderr*180*60*60*1000/math.pi
            self.raerrhms = self.raerrmas/(15*math.cos(self.decrad))

    def updateDec(self, newdecrad, newdecraderr=-1):
        self.decrad = newdecrad
        xxx, newdec = posradians2string(1, newdecrad)
        self.dec = newdec
        if newdecraderr > 0:
            self.decerrmas = newdecraderr*180*60*60*1000/math.pi

    def write(self, outputfile, overwrite=False):
        if os.path.exists(outputfile) and not overwrite:
            print("%s exists and overwrite is False - aborting" % outputfile)
            sys.exit()
        output = open(outputfile, "w")
        output.write("Pulsar %s: MJD %.4f: Frequency X, pol ?.:\n" % (self.object, self.mjd))
        output.write("%s%s\n" % ("Centre RA:".ljust(22), self.ra))
        output.write("%s%s\n" % ("Centre Dec:".ljust(22), self.dec))
        output.write("%s%.4f\n" % ("Flux (mJy):".ljust(22), self.flux))
        output.write("%s%.4f\n" % ("S/N:".ljust(22), self.snr))
        output.write("%s%.4f\n" % ("RA offset (mas)".ljust(22), 0.0))
        output.write("%s%.4f\n" % ("Dec offset (mas)".ljust(22), 0.0))
        output.write("%s%s\n" % ("Actual RA:".ljust(22), self.ra))
        output.write("%s%s\n" % ("Actual Dec:".ljust(22), self.dec))
        output.write("%s%.4fx%.4f at %.2f degrees\n" % ("Fit:".ljust(22), self.fitmin, self.fitmaj, self.fitpa))
        output.write("%s%.4f\n" % ("Est. RA error (mas):".ljust(22), self.raerrmas))
        output.write("%s%.4f\n" % ("Est. RA error (hms):".ljust(22), self.raerrhms))
        output.write("%s%.4f\n" % ("Est. Dec error (mas):".ljust(22), self.decerrmas))
        output.close()




#######################################################

# CALCULATING THE WEIGHTED MEAN POSITIONAL OFFSET

#######################################################


def weighted_avg_and_std(offsetval, weights):

    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.

    Parameters
    ----------
    offsetval : list
        List containing the offsets
    weights : list
        List containing the weights for the offsets

    Returns
    ----------
    lists
        wt_average : list containing the weighted average of the offsets
        std : list containing the weighted standard deviation of the offsets
    """

    # Define the number of points in the offset array
    npoints = len(offsetval)
    # Calculate the weighted average
    wt_average = np.average(offsetval, weights=weights ** 2)
    # Calculate the variance and standard deviation (fast and numerically precise method)
    variance = np.average((offsetval - wt_average) ** 2, weights=weights ** 2)
    # Standard deviation in an average-weighted point about the mean; i.e., the standard
    # deviation of the points about the mean
    std = np.sqrt(variance)
    # The error in the mean
    std_mean = std / np.sqrt(npoints - 1)

    return (wt_average, std_mean)


def plotter(
    ra_offsets,
    dec_offsets,
    ra_offset_unc,
    dec_offset_unc,
    wavg_ra_offset,
    wavg_dec_offset,
    wavg_ra_offset_unc,
    wavg_dec_offset_unc,
    catname,
):

    """
    Plot the RA, Dec offsets and the derived weighted mean offset for a given comparison

    Parameters
    ----------
    ra_offsets : list
        List containing the RA offsets
    dec_offsets : list
        List containing the Dec offsets
    ra_offset_unc : list
        List containing the RA offset uncertainties
    dec_offset_unc : list
        List containing the Dec offset uncertainties
    wavg_ra_offset : float
        Weighted average RA offset
    wavg_dec_offset : float
        Weighted average Dec offset
    wavg_ra_offset_unc : float
        Weighted average RA offset uncertainty
    wavg_dec_offset_unc : float
        Weighted average Dec offset uncertainty
    catname : str
        Name of the catalogue being used for comparison

    Returns
    ----------
    plot
        RA, Dec offsets for each source in the comparison
        Weighted average RA, Dec offset for all sources in comparison
    """

    print(
        f"RA weighted mean offset: {wavg_ra_offset:0.2f} +/- {wavg_ra_offset_unc:0.2f}"
    )
    print(
        f"Dec weighted mean offset: {wavg_dec_offset:0.2f} +/- {wavg_dec_offset_unc:0.2f}"
    )

    plt.figure(figsize=(11, 11))
    plt.title(
        fr"{args.frbtitletext} field source offsets (arcsec) from the {catname} catalogue positions: "
        + "\n"
        r"RA weighted mean offset: "
        + f"{wavg_ra_offset:0.2f}"
        + r" +/- "
        + f"{wavg_ra_offset_unc:0.2f}"
        + "\n"
        r"Dec weighted mean offset: "
        + f"{wavg_dec_offset:0.2f}"
        + r" +/- "
        + f"{wavg_dec_offset_unc:0.2f}"
    )

    for n in np.arange(len(ra_offsets)):
        plt.errorbar(
            ra_offsets[n],
            dec_offsets[n],
            yerr=dec_offset_unc[n],
            xerr=ra_offset_unc[n],
            fmt="k.",
            elinewidth=2,
            ecolor=col_bank[f"{n%16+1}"],
            capsize=2,
            label=askap_names[n],
        )

    plt.errorbar(
        wavg_ra_offset,
        wavg_dec_offset,
        xerr=wavg_ra_offset_unc,
        yerr=wavg_dec_offset_unc,
        fmt="kx",
        elinewidth=2.5,
        capsize=2,
        label="Weighted Mean Offset",
    )

    plt.xlabel("RA offset [arcsec]")
    plt.ylabel("Dec offset [arcsec]")
    plt.legend()
    plt.savefig(f"{args.frbtitletext}_field_offsets_from_{catname}.png")


#######################################################

# BOOTSTRAPPING

#######################################################


def bootstrap(data, samp_size, num_trials):

    boot_mean = []
    boot_68 = []

    n = 0
    while n < num_trials:

        samp = np.random.choice(data, size=samp_size, replace=True)

        samp_mean = np.mean(samp)
        boot_mean.append(samp_mean)

        n += 1

    mean_boot_mean = np.mean(boot_mean)
    # add histogram
    # bin histogram in 5% increments and then take bin with hightest value --> call this the most probable value
    # show median and mode and their +/- confidence intervals

    return mean_boot_mean


#######################################################

# CALCULATING THE CHI-SQUARED AND REDUCED CHI-SQUARED

#######################################################


def chi_sq(meas, model, err, num_srcs):

    """
    Return the chi-squared value for a given set of measurements and model.
    meas, model -- Numpy arrays with the same shape
    err, deg_free -- floats
    """

    deg_free = num_srcs - 1
    print("Number of degrees of freedom: ", deg_free)
    chi_squared = np.sum(((meas - model) / err) ** 2)
    chisq_red = chi_squared / deg_free

    return (chi_squared, chisq_red)


# ----------------------------------------------------------------------------------------------------------
#
# Calculate offsets along beam axes
#
# ----------------------------------------------------------------------------------------------------------

def rotated_offsets(catnameref, refras, refdecs, sigmarefras, sigmarefdecs, jmfitnamefile, frbfitsfile):
    ''' 
    This function calculates offsets long the the synthesized beam axes
    
    Written by AB based on AD's hardcoded script
    
    Inputs - 
           Name of reference catalogue
           List of reference RAs in arcsec
           List of reference DECs in arcsec
           Uncertainties on reference RAs in arcsec
           Uncertainties on reference DECs in arcsec
           Filename of list of relevant jmfit files
           Fits image containing the FRB - required to get the synthesized beam
           
    Returns -
           Chi-squared of fit
           DoF of fit
    
    Writes to disk - 
           Offsets along the beam axes and associated uncertainties
    '''

    #--- Reading BPA from the FRB FITS file

    frbfits  = fits.open(frbfitsfile)
    frbhdr   = frbfits[0].header
    frbbpa   = frbhdr['BPA']    
    print("BPA from %s is %.2f deg"%(frbfitsfile,frbbpa))
    posangle = frbbpa*np.pi/180
    frbfits.close()
    
    #--- Reading list of field sources

    reflines = open(jmfitnamefile).readlines()    
    nrefsrcs = len(reflines)
    print("Found %d sources"%(nrefsrcs))
    #print(refras,refdecs,sigmarefras,sigmarefdecs)

    #--- Initializing necessary lists and parameters

    xoffsets = []
    yoffsets = []
    totalxuncertainties = []
    totalyuncertainties = []
    dof = 0
    chisq = 0

    #---Looping over field sources

    for isrc in range(0,nrefsrcs):
        
        #--- Reading details of the field source 

        askapfile = reflines[isrc].split('\n')[0]
        askapresult = AstroObsResult("src"+str(isrc), "CRAFT", open(askapfile).readlines(), False)
        
        #--- Reading details of the same source from the reference catalogue
        
        raref = refras[isrc]/RAD_IN_ARCSEC
        decref = refdecs[isrc]/RAD_IN_ARCSEC
        sigmararef = sigmarefras[isrc]
        sigmadecref = sigmarefdecs[isrc]
        
        #--- Calculating offsets and associated uncertainties

        raraddiff = posdiff(raref, decref, askapresult.rarad, decref)
        if(askapresult.rarad > raref):
            raraddiff = -raraddiff
        
        decraddiff = posdiff(raref, decref, raref, askapresult.decrad)
        if(askapresult.decrad > decref):
            decraddiff = -decraddiff

        xdiff = raraddiff*np.sin(posangle) + decraddiff*np.cos(posangle)
        ydiff = raraddiff*np.cos(posangle) - decraddiff*np.sin(posangle)
        deltapa = posangle - np.pi*askapresult.fitpa / 180.0
        xuncertainty = np.sqrt((askapresult.fitmaj*np.cos(deltapa))**2 + (askapresult.fitmin*np.sin(deltapa))**2) / (2350 * askapresult.snr)
        yuncertainty = np.sqrt((askapresult.fitmaj*np.sin(deltapa))**2 + (askapresult.fitmin*np.cos(deltapa))**2) / (2350 * askapresult.snr)
            
        rauncertainty = askapresult.raerrmas/1000.
        decuncertainty = askapresult.decerrmas/1000.
        xoffsigma = (180*60*60/np.pi) * xdiff / np.sqrt(xuncertainty**2 + sigmararef**2) 
        # Just using RA uncertainty from reference, since beam is pretty circular
        yoffsigma = (180*60*60/np.pi) * ydiff / np.sqrt(yuncertainty**2 + sigmararef**2) 
        # Just using RA uncertainty from reference, since beam is pretty circular
        
        print (isrc, raraddiff*180*60*60/np.pi, decraddiff*180*60*60/np.pi, rauncertainty, decuncertainty, \
               xdiff*180*60*60/np.pi, ydiff*180*60*60/np.pi, xuncertainty, yuncertainty, \
               np.sqrt(rauncertainty**2 + decuncertainty**2)/np.sqrt(xuncertainty**2 + yuncertainty**2), \
               xoffsigma, yoffsigma)
        
        #--- Appending results to the list of offsets

        xoffsets.append(xdiff*180*60*60*1000/np.pi)
        yoffsets.append(ydiff*180*60*60*1000/np.pi)
        chisq += xoffsigma**2 + yoffsigma**2
        dof += 2
        totalxuncertainties.append(1000*np.sqrt(xuncertainty**2 + sigmararef**2))
        totalyuncertainties.append(1000*np.sqrt(yuncertainty**2 + sigmararef**2))
        
    print("For the case of zero offset, chi squared is", chisq, "for", dof, "degrees of freedom")
    
    #--- Writing the offsets to the  output file

    output = open("askap2"+catnameref+"_rotated_offsets.dat","w")
    
    for i in range(len(xoffsets)):
        output.write(str(xoffsets[i]) + " ")
    output.write("\n")

    for i in range(len(yoffsets)):
        output.write(str(yoffsets[i]) + " ")
    output.write("\n")

    for i in range(len(totalxuncertainties)):
        output.write(str(totalxuncertainties[i]) + " ")
    output.write("\n")

    for i in range(len(totalyuncertainties)):
        output.write(str(totalyuncertainties[i]) + " ")
    output.write("\n")
    output.close()  
    
    #--- Work done! Now return the chi-squared and the DoF in case it is useful
    
    return(chisq,dof)




#-----------------------------------------------------------------------------------

# The main script starts here

#####################################################################################

# DETERMINE OFFSETS AND UNCERTAINTIES

#####################################################################################

# Grab all the source info and put in a list

askap_srcinfo = grabsrcinfo(askap_posfile)

if args.first is not None:
    first_srcinfo = grabsrcinfo(first_posfile)
if args.nvss is not None:
    nvss_srcinfo = grabsrcinfo(nvss_posfile)
if args.racs is not None:
    racs_srcinfo = grabsrcinfo(racs_posfile)
if args.sumss is not None:
    sumss_srcinfo = grabsrcinfo(sumss_posfile)
if args.vlass is not None:
    vlass_srcinfo = grabsrcinfo(vlass_posfile)

# Extract the positions in arcsec and radians

askap_radec = getradec(askap_srcinfo)
askap_ra_arcsec = askap_radec["ra_arc"]
askap_dec_arcsec = askap_radec["dec_arc"]
askap_dec_rad = askap_radec["dec_rad"]

if args.first is not None:
    first_radec = getradec(first_srcinfo, first=True)
    first_ra_arcsec = first_radec["ra_arc"]
    first_dec_arcsec = first_radec["dec_arc"]
    first_dec_rad = first_radec["dec_rad"]
if args.nvss is not None:
    nvss_radec = getradec(nvss_srcinfo)
    nvss_ra_arcsec = nvss_radec["ra_arc"]
    nvss_dec_arcsec = nvss_radec["dec_arc"]
    nvss_dec_rad = nvss_radec["dec_rad"]
if args.racs is not None:
    racs_radec = getradec(racs_srcinfo)
    racs_ra_arcsec = racs_radec["ra_arc"]
    racs_dec_arcsec = racs_radec["dec_arc"]
    racs_dec_rad = racs_radec["dec_rad"]
if args.sumss is not None:
    sumss_radec = getradec(sumss_srcinfo)
    sumss_ra_arcsec = sumss_radec["ra_arc"]
    sumss_dec_arcsec = sumss_radec["dec_arc"]
    sumss_dec_rad = sumss_radec["dec_rad"]
if args.vlass is not None:
    vlass_radec = getradec(vlass_srcinfo)
    vlass_ra_arcsec = vlass_radec["ra_arc"]
    vlass_dec_arcsec = vlass_radec["dec_arc"]
    vlass_dec_rad = vlass_radec["dec_rad"]


# Get uncertainties for RA and Dec in arcsec

askap_radec_unc = getradec_unc(askap_srcinfo, askap_dec_rad, askap=True)
askap_ra_unc_arcsec = askap_radec_unc["ra_unc_arcsec"]
askap_dec_unc_arcsec = askap_radec_unc["dec_unc_arcsec"]

if args.first is not None:
    first_radec_unc = getradec_unc(first_srcinfo, first_dec_rad, first=True)
    first_ra_unc_arcsec = first_radec_unc["ra_unc_arcsec"]
    first_dec_unc_arcsec = first_radec_unc["dec_unc_arcsec"]
if args.nvss is not None:
    nvss_radec_unc = getradec_unc(nvss_srcinfo, nvss_dec_rad, nvss=True)
    nvss_ra_unc_arcsec = nvss_radec_unc["ra_unc_arcsec"]
    nvss_dec_unc_arcsec = nvss_radec_unc["dec_unc_arcsec"]
if args.racs is not None:
    racs_radec_unc = getradec_unc(racs_srcinfo, racs_dec_rad, racs=True)
    racs_ra_unc_arcsec = racs_radec_unc["ra_unc_arcsec"]
    racs_dec_unc_arcsec = racs_radec_unc["dec_unc_arcsec"]
if args.sumss is not None:
    sumss_radec_unc = getradec_unc(sumss_srcinfo, sumss_dec_rad, sumss=True)
    sumss_ra_unc_arcsec = sumss_radec_unc["ra_unc_arcsec"]
    sumss_dec_unc_arcsec = sumss_radec_unc["dec_unc_arcsec"]
if args.vlass is not None:
    vlass_radec_unc = getradec_unc(vlass_srcinfo, vlass_dec_rad, vlass=True)
    vlass_ra_unc_arcsec = vlass_radec_unc["ra_unc_arcsec"]
    vlass_dec_unc_arcsec = vlass_radec_unc["dec_unc_arcsec"]

# Determine the offsets and their total uncertainties between the ASKAP and catalogue source positions,
# and save the offsets and uncertainties to a .dat file for use with weighted_multi_image_fit.py

if args.first is not None:
    # Offsets
    askap2first_offsets = getoffsets(
        askap_ra_arcsec,
        askap_dec_arcsec,
        first_ra_arcsec,
        first_dec_arcsec,
        first_dec_rad,
    )
    askap2first_offsets_ra = askap2first_offsets["ra_offset_arcsec"]
    askap2first_offsets_dec = askap2first_offsets["dec_offset_arcsec"]
    # Total uncertainties
    askap2first_offsets_unc = gettotaloffsetunc(
        askap_ra_unc_arcsec,
        askap_dec_unc_arcsec,
        first_ra_unc_arcsec,
        first_dec_unc_arcsec,
    )
    askap2first_offsets_unc_ra = askap2first_offsets_unc["totra_unc"]
    askap2first_offsets_unc_dec = askap2first_offsets_unc["totdec_unc"]
    # Format into strings
    offsets_ra_first_str = " ".join(
        [
            str(askap2first_offsets_ra[i])
            for i in np.arange(len(askap2first_offsets_ra))
        ]
    )
    offsets_dec_first_str = " ".join(
        [
            str(askap2first_offsets_dec[i])
            for i in np.arange(len(askap2first_offsets_dec))
        ]
    )
    offsets_ra_unc_first_str = " ".join(
        [
            str(askap2first_offsets_unc_ra[i])
            for i in np.arange(len(askap2first_offsets_unc_ra))
        ]
    )
    offsets_dec_unc_first_str = " ".join(
        [
            str(askap2first_offsets_unc_dec[i])
            for i in np.arange(len(askap2first_offsets_unc_dec))
        ]
    )
    # Save to .dat file
    offset_unc_out = open("askap2first_offsets_unc.dat", "w")
    offset_unc_out.write(f"{offsets_ra_first_str}\n")
    offset_unc_out.write(f"{offsets_dec_first_str}\n")
    offset_unc_out.write(f"{offsets_ra_unc_first_str}\n")
    offset_unc_out.write(f"{offsets_dec_unc_first_str}\n")
    offset_unc_out.close()

if args.nvss is not None:
    # Offsets
    askap2nvss_offsets = getoffsets(
        askap_ra_arcsec,
        askap_dec_arcsec,
        nvss_ra_arcsec,
        nvss_dec_arcsec,
        nvss_dec_rad,
    )
    askap2nvss_offsets_ra = askap2nvss_offsets["ra_offset_arcsec"]
    askap2nvss_offsets_dec = askap2nvss_offsets["dec_offset_arcsec"]
    # Total uncertainties
    askap2nvss_offsets_unc = gettotaloffsetunc(
        askap_ra_unc_arcsec,
        askap_dec_unc_arcsec,
        nvss_ra_unc_arcsec,
        nvss_dec_unc_arcsec,
    )
    askap2nvss_offsets_unc_ra = askap2nvss_offsets_unc["totra_unc"]
    askap2nvss_offsets_unc_dec = askap2nvss_offsets_unc["totdec_unc"]
    # Format into strings
    offsets_ra_nvss_str = " ".join(
        [
            str(askap2nvss_offsets_ra[i])
            for i in np.arange(len(askap2nvss_offsets_ra))
        ]
    )
    offsets_dec_nvss_str = " ".join(
        [
            str(askap2nvss_offsets_dec[i])
            for i in np.arange(len(askap2nvss_offsets_dec))
        ]
    )
    offsets_ra_unc_nvss_str = " ".join(
        [
            str(askap2nvss_offsets_unc_ra[i])
            for i in np.arange(len(askap2nvss_offsets_unc_ra))
        ]
    )
    offsets_dec_unc_nvss_str = " ".join(
        [
            str(askap2nvss_offsets_unc_dec[i])
            for i in np.arange(len(askap2nvss_offsets_unc_dec))
        ]
    )
    # Save to .dat file
    offset_unc_out = open("askap2nvss_offsets_unc.dat", "w")
    offset_unc_out.write(f"{offsets_ra_nvss_str}\n")
    offset_unc_out.write(f"{offsets_dec_nvss_str}\n")
    offset_unc_out.write(f"{offsets_ra_unc_nvss_str}\n")
    offset_unc_out.write(f"{offsets_dec_unc_nvss_str}\n")
    offset_unc_out.close()

if args.racs is not None:
    # Offsets
    askap2racs_offsets = getoffsets(
        askap_ra_arcsec,
        askap_dec_arcsec,
        racs_ra_arcsec,
        racs_dec_arcsec,
        racs_dec_rad,
    )
    askap2racs_offsets_ra = askap2racs_offsets["ra_offset_arcsec"]
    askap2racs_offsets_dec = askap2racs_offsets["dec_offset_arcsec"]
    # Total uncertainties
    askap2racs_offsets_unc = gettotaloffsetunc(
        askap_ra_unc_arcsec,
        askap_dec_unc_arcsec,
        racs_ra_unc_arcsec,
        racs_dec_unc_arcsec,
    )
    askap2racs_offsets_unc_ra = askap2racs_offsets_unc["totra_unc"]
    askap2racs_offsets_unc_dec = askap2racs_offsets_unc["totdec_unc"]
    # Format into strings
    offsets_ra_racs_str = " ".join(
        [
            str(askap2racs_offsets_ra[i])
            for i in np.arange(len(askap2racs_offsets_ra))
        ]
    )
    offsets_dec_racs_str = " ".join(
        [
            str(askap2racs_offsets_dec[i])
            for i in np.arange(len(askap2racs_offsets_dec))
        ]
    )
    offsets_ra_unc_racs_str = " ".join(
        [
            str(askap2racs_offsets_unc_ra[i])
            for i in np.arange(len(askap2racs_offsets_unc_ra))
        ]
    )
    offsets_dec_unc_racs_str = " ".join(
        [
            str(askap2racs_offsets_unc_dec[i])
            for i in np.arange(len(askap2racs_offsets_unc_dec))
        ]
    )
    # Save to .dat file
    offset_unc_out = open("askap2racs_offsets_unc.dat", "w")
    offset_unc_out.write(f"{offsets_ra_racs_str}\n")
    offset_unc_out.write(f"{offsets_dec_racs_str}\n")
    offset_unc_out.write(f"{offsets_ra_unc_racs_str}\n")
    offset_unc_out.write(f"{offsets_dec_unc_racs_str}\n")
    offset_unc_out.close()

if args.sumss is not None:
    # Offsets
    askap2sumss_offsets = getoffsets(
        askap_ra_arcsec,
        askap_dec_arcsec,
        sumss_ra_arcsec,
        sumss_dec_arcsec,
        sumss_dec_rad,
    )
    askap2sumss_offsets_ra = askap2sumss_offsets["ra_offset_arcsec"]
    askap2sumss_offsets_dec = askap2sumss_offsets["dec_offset_arcsec"]
    # Total uncertainties
    askap2sumss_offsets_unc = gettotaloffsetunc(
        askap_ra_unc_arcsec,
        askap_dec_unc_arcsec,
        sumss_ra_unc_arcsec,
        sumss_dec_unc_arcsec,
    )
    askap2sumss_offsets_unc_ra = askap2sumss_offsets_unc["totra_unc"]
    askap2sumss_offsets_unc_dec = askap2sumss_offsets_unc["totdec_unc"]
    # Format into strings
    offsets_ra_sumss_str = " ".join(
        [
            str(askap2sumss_offsets_ra[i])
            for i in np.arange(len(askap2sumss_offsets_ra))
        ]
    )
    offsets_dec_sumss_str = " ".join(
        [
            str(askap2sumss_offsets_dec[i])
            for i in np.arange(len(askap2sumss_offsets_dec))
        ]
    )
    offsets_ra_unc_sumss_str = " ".join(
        [
            str(askap2sumss_offsets_unc_ra[i])
            for i in np.arange(len(askap2sumss_offsets_unc_ra))
        ]
    )
    offsets_dec_unc_sumss_str = " ".join(
        [
            str(askap2sumss_offsets_unc_dec[i])
            for i in np.arange(len(askap2sumss_offsets_unc_dec))
        ]
    )
    # Save to .dat file
    offset_unc_out = open("askap2sumss_offsets_unc.dat", "w")
    offset_unc_out.write(f"{offsets_ra_sumss_str}\n")
    offset_unc_out.write(f"{offsets_dec_sumss_str}\n")
    offset_unc_out.write(f"{offsets_ra_unc_sumss_str}\n")
    offset_unc_out.write(f"{offsets_dec_unc_sumss_str}\n")
    offset_unc_out.close()

if args.vlass is not None:
    # Offsets
    askap2vlass_offsets = getoffsets(
        askap_ra_arcsec,
        askap_dec_arcsec,
        vlass_ra_arcsec,
        vlass_dec_arcsec,
        vlass_dec_rad,
    )
    askap2vlass_offsets_ra = askap2vlass_offsets["ra_offset_arcsec"]
    askap2vlass_offsets_dec = askap2vlass_offsets["dec_offset_arcsec"]
    # Total uncertainties
    askap2vlass_offsets_unc = gettotaloffsetunc(
        askap_ra_unc_arcsec,
        askap_dec_unc_arcsec,
        vlass_ra_unc_arcsec,
        vlass_dec_unc_arcsec,
    )
    askap2vlass_offsets_unc_ra = askap2vlass_offsets_unc["totra_unc"]
    askap2vlass_offsets_unc_dec = askap2vlass_offsets_unc["totdec_unc"]
    # Format into strings
    offsets_ra_vlass_str = " ".join(
        [
            str(askap2vlass_offsets_ra[i])
            for i in np.arange(len(askap2vlass_offsets_ra))
        ]
    )
    offsets_dec_vlass_str = " ".join(
        [
            str(askap2vlass_offsets_dec[i])
            for i in np.arange(len(askap2vlass_offsets_dec))
        ]
    )
    offsets_ra_unc_vlass_str = " ".join(
        [
            str(askap2vlass_offsets_unc_ra[i])
            for i in np.arange(len(askap2vlass_offsets_unc_ra))
        ]
    )
    offsets_dec_unc_vlass_str = " ".join(
        [
            str(askap2vlass_offsets_unc_dec[i])
            for i in np.arange(len(askap2vlass_offsets_unc_dec))
        ]
    )
    # Save to .dat file
    offset_unc_out = open("askap2vlass_offsets_unc.dat", "w")
    offset_unc_out.write(f"{offsets_ra_vlass_str}\n")
    offset_unc_out.write(f"{offsets_dec_vlass_str}\n")
    offset_unc_out.write(f"{offsets_ra_unc_vlass_str}\n")
    offset_unc_out.write(f"{offsets_dec_unc_vlass_str}\n")
    offset_unc_out.close()


# Weight the RA and Dec offsets by their total errors
if args.first is not None:
    raweight_first = 1.0 / np.array(askap2first_offsets_unc_ra)
    decweight_first = 1.0 / np.array(askap2first_offsets_unc_dec)
if args.nvss is not None:
    raweight_nvss = 1.0 / np.array(askap2nvss_offsets_unc_ra)
    decweight_nvss = 1.0 / np.array(askap2nvss_offsets_unc_dec)
if args.racs is not None:
    raweight_racs = 1.0 / np.array(askap2racs_offsets_unc_ra)
    decweight_racs = 1.0 / np.array(askap2racs_offsets_unc_dec)
if args.sumss is not None:
    raweight_sumss = 1.0 / np.array(askap2sumss_offsets_unc_ra)
    decweight_sumss = 1.0 / np.array(askap2sumss_offsets_unc_dec)
if args.vlass is not None:
    raweight_vlass = 1.0 / np.array(askap2vlass_offsets_unc_ra)
    decweight_vlass = 1.0 / np.array(askap2vlass_offsets_unc_dec)

# Calculate the weighted mean offset in RA and Dec
if args.first is not None:
    wmean_raoffset_first, wmean_raoffsetunc_first = weighted_avg_and_std(
        askap2first_offsets_ra, raweight_first
    )
    wmean_decoffset_first, wmean_decoffsetunc_first = weighted_avg_and_std(
        askap2first_offsets_dec, decweight_first
    )
if args.nvss is not None:
    wmean_raoffset_nvss, wmean_raoffsetunc_nvss = weighted_avg_and_std(
        askap2nvss_offsets_ra, raweight_nvss
    )
    wmean_decoffset_nvss, wmean_decoffsetunc_nvss = weighted_avg_and_std(
        askap2nvss_offsets_dec, decweight_nvss
    )
if args.racs is not None:
    wmean_raoffset_racs, wmean_raoffsetunc_racs = weighted_avg_and_std(
        askap2racs_offsets_ra, raweight_racs
    )
    wmean_decoffset_racs, wmean_decoffsetunc_racs = weighted_avg_and_std(
        askap2racs_offsets_dec, decweight_racs
    )
if args.sumss is not None:
    wmean_raoffset_sumss, wmean_raoffsetunc_sumss = weighted_avg_and_std(
        askap2sumss_offsets_ra, raweight_sumss
    )
    wmean_decoffset_sumss, wmean_decoffsetunc_sumss = weighted_avg_and_std(
        askap2sumss_offsets_dec, decweight_sumss
    )
if args.vlass is not None:
    wmean_raoffset_vlass, wmean_raoffsetunc_vlass = weighted_avg_and_std(
        askap2vlass_offsets_ra, raweight_vlass
    )
    wmean_decoffset_vlass, wmean_decoffsetunc_vlass = weighted_avg_and_std(
        askap2vlass_offsets_dec, decweight_vlass
    )

# Grab ASKAP source names
askap_names = grabsrcinfo(args.askapnames)
askap_names = [askap_names[s].strip("\n") for s in np.arange(len(askap_names))]


# Plot the comparisons and save the plots
if args.first is not None:
    plotter(
        askap2first_offsets_ra,
        askap2first_offsets_dec,
        askap2first_offsets_unc_ra,
        askap2first_offsets_unc_dec,
        wmean_raoffset_first,
        wmean_decoffset_first,
        wmean_raoffsetunc_first,
        wmean_decoffsetunc_first,
        catname="FIRST",
    )
if args.nvss is not None:
    plotter(
        askap2nvss_offsets_ra,
        askap2nvss_offsets_dec,
        askap2nvss_offsets_unc_ra,
        askap2nvss_offsets_unc_dec,
        wmean_raoffset_nvss,
        wmean_decoffset_nvss,
        wmean_raoffsetunc_nvss,
        wmean_decoffsetunc_nvss,
        catname="NVSS",
    )
if args.racs is not None:
    plotter(
        askap2racs_offsets_ra,
        askap2racs_offsets_dec,
        askap2racs_offsets_unc_ra,
        askap2racs_offsets_unc_dec,
        wmean_raoffset_racs,
        wmean_decoffset_racs,
        wmean_raoffsetunc_racs,
        wmean_decoffsetunc_racs,
        catname="racs",
    )
if args.sumss is not None:
    plotter(
        askap2sumss_offsets_ra,
        askap2sumss_offsets_dec,
        askap2sumss_offsets_unc_ra,
        askap2sumss_offsets_unc_dec,
        wmean_raoffset_sumss,
        wmean_decoffset_sumss,
        wmean_raoffsetunc_sumss,
        wmean_decoffsetunc_sumss,
        catname="SUMSS",
    )
if args.vlass is not None:
    plotter(
        askap2vlass_offsets_ra,
        askap2vlass_offsets_dec,
        askap2vlass_offsets_unc_ra,
        askap2vlass_offsets_unc_dec,
        wmean_raoffset_vlass,
        wmean_decoffset_vlass,
        wmean_raoffsetunc_vlass,
        wmean_decoffsetunc_vlass,
        catname="VLASS",
    )

# Calculate the numbers of 'points' of measurement; i.e., the number of sources
# num_sources = len(askapnames)

# Two test models:
# Test 1: A zero offset in position
# pos_model_zero_ra = np.zeros(np.shape(ra_offset))
# pos_model_zero_dec = np.zeros(np.shape(ra_offset))
# chi_sq_zero_ra, chisq_red_zero_ra = chi_sq(ra_offset, pos_model_zero_ra, totra_err, num_sources)
# chi_sq_zero_dec, chisq_red_zero_dec = chi_sq(dec_offset, pos_model_zero_dec, totdec_err, num_sources)
# Test 2: Offsets given by the above weighted mean offsets
# pos_model_wm_ra = np.full(np.shape(ra_offset),wmean_raoffset)
# pos_model_wm_dec = np.full(np.shape(dec_offset),wmean_decoffset)
# chi_sq_wm_ra, chisq_red_wm_ra = chi_sq(ra_offset, pos_model_wm_ra, totra_err, num_sources)
# chi_sq_wm_dec, chisq_red_wm_dec = chi_sq(dec_offset, pos_model_wm_dec, totdec_err, num_sources)

# print("chi-squared, reduced chi-squared RA (zero offset model): ", chi_sq_zero_ra, chisq_red_zero_ra)
# print("chi-squared, reduced chi-squared Dec (zero offset model): ", chi_sq_zero_dec, chisq_red_zero_dec)
# print("chi-squared, reduced chi-squared RA (weighted mean offset model): ", chi_sq_wm_ra, chisq_red_wm_ra)
# print("chi-squared, reduced chi-squared Dec (weighted mean offset model): ", chi_sq_wm_dec, chisq_red_wm_dec)
# Calculate the numbers of 'points' of measurement; i.e., the number of sources
# num_sources = len(askapnames)

# Two test models:
# Test 1: A zero offset in position
# pos_model_zero_ra = np.zeros(np.shape(ra_offset))
# pos_model_zero_dec = np.zeros(np.shape(ra_offset))
# chi_sq_zero_ra, chisq_red_zero_ra = chi_sq(ra_offset, pos_model_zero_ra, totra_err, num_sources)
# chi_sq_zero_dec, chisq_red_zero_dec = chi_sq(dec_offset, pos_model_zero_dec, totdec_err, num_sources)
# Test 2: Offsets given by the above weighted mean offsets
# pos_model_wm_ra = np.full(np.shape(ra_offset),wmean_raoffset)
# pos_model_wm_dec = np.full(np.shape(dec_offset),wmean_decoffset)
# chi_sq_wm_ra, chisq_red_wm_ra = chi_sq(ra_offset, pos_model_wm_ra, totra_err, num_sources)
# chi_sq_wm_dec, chisq_red_wm_dec = chi_sq(dec_offset, pos_model_wm_dec, totdec_err, num_sources)

# print("chi-squared, reduced chi-squared RA (zero offset model): ", chi_sq_zero_ra, chisq_red_zero_ra)
# print("chi-squared, reduced chi-squared Dec (zero offset model): ", chi_sq_zero_dec, chisq_red_zero_dec)
# print("chi-squared, reduced chi-squared RA (weighted mean offset model): ", chi_sq_wm_ra, chisq_red_wm_ra)
# print("chi-squared, reduced chi-squared Dec (weighted mean offset model): ", chi_sq_wm_dec, chisq_red_wm_dec)

#--------------------------------------------------------------------------------------------------------------
# Alright, now calculate offsets along the beam axes

if args.first is not None:
    dummyres = rotated_offsets('first',first_ra_arcsec,first_dec_arcsec,first_ra_unc_arcsec,first_dec_unc_arcsec,args.jmfitnames,args.fieldfits)
    print(dummyres)

if args.nvss is not None:
    dummyres = rotated_offsets('nvss',nvss_ra_arcsec,nvss_dec_arcsec,nvss_ra_unc_arcsec,nvss_dec_unc_arcsec,args.jmfitnames,args.fieldfits)
    print(dummyres)

if args.racs is not None:
    dummyres = rotated_offsets('racs',racs_ra_arcsec,racs_dec_arcsec,racs_ra_unc_arcsec,racs_dec_unc_arcsec,args.jmfitnames,args.fieldfits)
    print(dummyres)

if args.sumss is not None:
    dummyres = rotated_offsets('sumss',sumss_ra_arcsec,sumss_dec_arcsec,sumss_ra_unc_arcsec,sumss_dec_unc_arcsec,args.jmfitnames,args.fieldfits)
    print(dummyres)

if args.vlass is not None:
    dummyres = rotated_offsets('vlass',racs_ra_arcsec,racs_dec_arcsec,racs_ra_unc_arcsec,racs_dec_unc_arcsec,args.jmfitnames,args.fieldfits)
    print(dummyres)
