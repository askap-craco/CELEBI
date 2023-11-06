###############################################################################
# craftpy2 stage 2.6: apply offset
# Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
# Modified by AB based on AD's hardcoded script: 2023-03-16
#
# Applies the systematic offset calculated from field source positions and
# applies it to the FRB position
###############################################################################

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import os
import math
import numpy as np
from astropy import units as un
from astropy.io import fits
from astropy.coordinates import SkyCoord as sc

parser = ArgumentParser(
    description="Apply the systematic offset calculated "
    "from field source positions and apply it "
    "to the FRB position.",
    formatter_class=ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--frbname", help="FRB name string")
parser.add_argument("--frb", help="Pre-offset FRB JMFIT stats file")
parser.add_argument("--offset", help="Offset file")
parser.add_argument("--doffset", help="Detailed offset file")
parser.add_argument("--frbfits", help="FITS image containing the FRB")
parser.add_argument("--hpfits", help="Healpix FITS map for FRB position")

args = parser.parse_args()

SYSTEMATIC_SCALE_FACTOR = 1.79
RADtoDEG = 180.0/np.pi

#--- Reading BPA from the FRB FITS file

frbfits  = fits.open(args.frbfits)
frbhdr   = frbfits[0].header
frbbpa   = frbhdr['BPA']    
#print("BPA from %s is %.2f deg"%(args.frbfits,frbbpa))
posangle = frbbpa*np.pi/180
frbfits.close()

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


class Coord:
    def __init__(self, stats):
        fields = {}
        f = open(stats)
        for line in f:
            line = line[:-1].split(": ")  # [:-1] to trim newline
            fields[line[0]] = line[1].strip()  # strip extra whitespace
        f.close()

        self.ra_hms = fields["Actual RA"]  # hms
        self.ra_err = float(fields["Est. RA error (mas)"]) / 1e3  # arcseconds
        self.dec_dms = fields["Actual Dec"]  # dms
        self.dec_err = (
            float(fields["Est. Dec error (mas)"]) / 1e3
        )  # arcseconds


FRB = Coord(args.frb)
offset_ra, offset_ra_err, offset_dec, offset_dec_err = np.loadtxt(args.offset)

frb_pos = sc(FRB.ra_hms, FRB.dec_dms, unit=("hourangle", "deg"))

final_pos = sc(
    frb_pos.ra + offset_ra*un.arcsecond/np.cos(frb_pos.dec.to(un.rad)),
    frb_pos.dec + offset_dec*un.arcsecond
)

final_err_ra = (FRB.ra_err ** 2 + offset_ra_err ** 2) ** 0.5
final_err_dec = (FRB.dec_err ** 2 + offset_dec_err ** 2) ** 0.5


print("FINAL FRB POSITION (in the Old System)")
print(
    f'RA:\t{final_pos.ra.to_string("hour")}'
    f'\t+-{final_err_ra}'
)
print(
    f'Dec:\t{final_pos.dec.to_string("deg")}'
    f'\t+-{final_err_dec}'
)

# Now correct for offsets along the beam axes

frbresult       = AstroObsResult("FRB", "CRAFT", open(args.frb).readlines(), False)
lines           = open(args.doffset).readlines()

resultline = ""
for i in range(len(lines)):
    if "########## H1: offsets" in lines[i]:
        resultline = lines[i+3]
resline1  = resultline.replace("("," ")
resline2  = resline1.replace(")"," ")
resline3  = resline2.replace("["," ")
resline4  = resline3.replace("]"," ")
splitline = resline4.split()

majoffset = float(splitline[1])
majsigma  = float(splitline[4])*SYSTEMATIC_SCALE_FACTOR
minoffset = float(splitline[5])
minsigma  = float(splitline[8])*SYSTEMATIC_SCALE_FACTOR

raoffset  = majoffset*np.sin(posangle) + minoffset*np.cos(posangle)
decoffset = majoffset*np.cos(posangle) - minoffset*np.sin(posangle)
frbresult.updateRA(frbresult.rarad + raoffset*np.pi/(180*60*60*1000))
frbresult.updateDec(frbresult.decrad + decoffset*np.pi/(180*60*60*1000))
deltapa = posangle - np.pi*frbresult.fitpa / 180.0
frbstatmaj_unc = np.sqrt((frbresult.fitmaj*np.cos(deltapa))**2 + (frbresult.fitmin*np.sin(deltapa))**2) / (2.350 * frbresult.snr)
frbstatmin_unc = np.sqrt((frbresult.fitmaj*np.sin(deltapa))**2 + (frbresult.fitmin*np.cos(deltapa))**2) / (2.350 * frbresult.snr)

totalmaj_unc = np.sqrt(frbstatmaj_unc**2 + majsigma**2)
totalmin_unc = np.sqrt(frbstatmin_unc**2 + minsigma**2)

raerror   = np.sqrt((totalmaj_unc*np.sin(posangle))**2 + (totalmin_unc*np.cos(posangle))**2)
decerror  = np.sqrt((totalmaj_unc*np.cos(posangle))**2 + (totalmin_unc*np.sin(posangle))**2)

print("\n*********************************************\n")
print("Offsets along the beam direction\n")
print("Offsets       -- MAJOR = %.2f arcsec, MINOR = %.2f arcsec, PA = %.2f deg"%(majoffset/1000.0, minoffset/1000.0, posangle*180.0/np.pi))
print("Uncertainties -- MAJOR = %.2f arcsec, MINOR = %.2f arcsec"%(majsigma/1000.0,minsigma/1000.0))
print("\nRA, Dec offset is {0:.3f} arcsec, {1:.3f} arcsec".format(raoffset/1000.0, decoffset/1000.0))

print("\n*********************************************\n")
print("FINAL FRB POSITION for offsets along beam direction\n")

print("RA     {0} +/- {1:.3f}".format(frbresult.ra,raerror/1000.0))
print("DEC    {0} +/- {1:.3f}\n".format(frbresult.dec,decerror/1000.0))
print("Uncertainty ellipse\n")
print("Major axis = {0:.3f} arcsec, Minor axis = {1:.3f} arcsec, PA = {2:.2f} deg\n".format(totalmaj_unc/1000.0, totalmin_unc/1000.0, posangle*180/np.pi)) 
print("The systematic uncertainty (which was {0:.2f} x {1:.2f}) had a scale factor of {2} applied".format(majsigma/1e3, minsigma/1e3, SYSTEMATIC_SCALE_FACTOR))

# Comparing errors for debuggin purpose
#print("\n\nComparing uncertainties")
#print("\nWITHOUT rotated offsets\n")
#print("FRB RA error = %.5f   FRB Dec error = %.5f    Offset RA error = %.5f   Offset Dec error = %.5f"%(FRB.ra_err,FRB.dec_err,offset_ra_err,offset_dec_err))
#print("\nWITH rotated offsets\n")
#print("FRB MAJ error= %.5f   FRB MIN error = %.5f    Offset MAJ error= %.5f   Offset MIN error = %.5f"%(frbstatmaj_unc/1e3,frbstatmin_unc/1e3,majsigma/1e3,minsigma/1e3))

# Now generate a Healpix map that Prof X wants
print("\nGenerating Healpix map using craco_fu_hp .....")
#print("export PATH=$PATH:/fred/oz002/askap/craft/craco/anaconda/anaconda3/bin; unset PYTHONPATH; conda run -n cracofubase craco_fu_hp FRB{0} --coord {1},{2} --siga {3} --sigb {4} --PA {5} --clobber --outfile {6}".format(args.frbname,frbresult.rarad*RADtoDEG,frbresult.decrad*RADtoDEG,totalmaj_unc/1000.0,totalmin_unc/1000.0,posangle*180.0/np.pi,args.hpfits)) 

os.system("export PATH=/fred/oz002/askap/craft/craco/anaconda/anaconda3/bin:$PATH; unset PYTHONPATH; conda run -n cracofubase craco_fu_hp FRB{0} --coord {1},{2} --siga {3} --sigb {4} --PA {5} --clobber --outfile {6}".format(args.frbname,frbresult.rarad*RADtoDEG,frbresult.decrad*RADtoDEG,totalmaj_unc/1000.0,totalmin_unc/1000.0,posangle*180.0/np.pi,args.hpfits)) 




