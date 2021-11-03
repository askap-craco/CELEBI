import glob
import sys

data_path = sys.argv[1]
hdrfiles = glob.glob(f"{data_path}/ak*/*/*hdr")

# initialise to ensure scope
samprate = None
nsamps = None
thisMJD = None
startmjd = 10000000

# Heavily inspired by difx's vcraft2obs.py
for filename in hdrfiles:
    with open(filename) as file:
        for line in file:
            if line.split()[0] == "SAMP_RATE":
                samprate = float(line.split()[1])
            if line.split()[0] == "NSAMPS_REQUEST":
                nsamps = int(line.split()[1])
            if line.split()[0] == "TRIGGER_MJD":
                thisMJD = float(line.split()[1])

    new_startmjd = thisMJD - nsamps / (samprate * 86400)
    startmjd = new_startmjd if new_startmjd < startmjd else startmjd

print(startmjd)
