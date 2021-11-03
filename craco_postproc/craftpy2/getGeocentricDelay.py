#!/usr/bin/env python3

import glob
import math
import os
import sys

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

askap_lat = -26.697
askap_lon = 116.631
askap_height = 361
askap_radius = 6374217  # distance from geocentre
c = 299792458.0

if not len(sys.argv) == 3:
    print("Usage: %s <rawdata directory> <snoopy file>" % sys.argv[0])
    print("Raw data directory should contain akXX/beamXX/*.vcraft.hdr files")
    sys.exit()


rawdata = sys.argv[1]
snoopy_file = sys.argv[2]

if not os.path.exists(snoopy_file):
    print((snoopy_file, "doesn't exist"))
    sys.exit()


# Parse the snoopy file
with open(snoopy_file) as snoopyin:
    snoopylines = snoopyin.readlines()

nocommentlines = []
for line in snoopylines:
    print(line)
    if len(line) > 1 and not line[0] == "#":
        nocommentlines.append(line)
        print(("Snoopy info", nocommentlines))
if len(nocommentlines) != 1:
    print("ERROR: No information found")
    sys.exit()

splitline = nocommentlines[0].split()
mjd = float(splitline[7])


# Grab all the headers and find the lowest frequency and the rough FRB direction
antennas = sorted(glob.glob(sys.argv[1] + "/ak??/"))
if len(antennas) == 0:
    print("No antennas found!")
    sys.exit()
beams = glob.glob(antennas[0] + "/beam??/")
if len(beams) == 0:
    print(("No beams found in", antennas[0], "(in /beam??/)"))
    sys.exit()
hdrfiles = glob.glob(beams[0] + "/*hdr")
if len(hdrfiles) == 0:
    print("No vcraft header files found!")
    sys.exit()
lowestfreq = 99999999
for hf in hdrfiles:
    with open(hf) as headerin:
        lines = headerin.readlines()
        for line in lines:
            if "FREQ" in line:
                freqs = line.split()[1].split(",")
                for f in freqs:
                    if float(f) < lowestfreq:
                        lowestfreq = float(f)
            if "BEAM_RA" in line:
                beamra = float(line.split()[1])
            if "BEAM_DEC" in line:
                beamdec = float(line.split()[1])
            if "TRIGGER_MJD" in line:
                triggermjd = float(line.split()[1])
            if "SAMP_RATE" in line:
                samprate = float(line.split()[1])
            if "NSAMPS_REQUEST" in line:
                nsamps = int(line.split()[1])

            # if "ANT_EL" in line:
            #    ant_el_rad = float(line.split()[1])*math.pi/180.0

frb = SkyCoord(beamra, beamdec, unit="deg")
askap = EarthLocation(
    lon=askap_lon * u.deg, lat=askap_lat * u.deg, height=askap_height * u.m
)
t = Time(triggermjd, format="mjd")
altaz = frb.transform_to(AltAz(obstime=t, location=askap))
geocentricdelay = math.sin(altaz.alt.radian) * askap_radius / c

# Would be much simpler to just get the elevation from the vcraft file!!
# But this is for the centre of the PAF, so could be a bit off vs the elevation of the actual beam
# (which is itself already a bit off vs the actual FRB, but within half a degree which is good enough)
# geocentricdelay = math.sin(ant_el_rad) * askap_radius / c

# Now figure out what the correlation start time will be:
corrstartmjd = math.floor(triggermjd * 86400 - nsamps / samprate) / 86400.0

print(("Geocentric delay is", geocentricdelay))
print(("Lowest frequency is", lowestfreq))

print(
    f"snoopylog2frbgate.py -f {lowestfreq:.3f} --timediff {geocentricdelay * 1e3:.3f} --corrstartmjd {corrstartmjd:.9f} {snoopy_file}"
)
