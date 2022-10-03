# Write a profile as a function of MJD

import glob
import sys
import numpy as np

# args: data_path I_50us_dynspec 50us_crop_start_s.txt

data_path = sys.argv[1]
hdrfiles = glob.glob(f'{data_path}/ak*/*/*hdr')

# STEP 1: find data start time == last start MJD in headers
# initialise to ensure scope
samprate = None
nsamps = None
thisMJD = None
lastMJD = 0

# Heavily inspired by difx's vcraft2obs.py
for filename in hdrfiles:
    with open(filename) as file:
        for line in file:
            if line.split()[0] == 'SAMP_RATE':
                samprate = float(line.split()[1])
            if line.split()[0] == 'NSAMPS_REQUEST':
                nsamps = int(line.split()[1])
            if line.split()[0] == 'TRIGGER_MJD':
                thisMJD = float(line.split()[1])
            
    new_lastMJD = thisMJD - nsamps/(samprate*86400)
    lastMJD = new_lastMJD if new_lastMJD > lastMJD else lastMJD


# STEP 2: fscrunch dynamic spectrum to get profile
I = np.sum(np.load(sys.argv[2]), axis=0)

# STEP 3: create time axis, assuming 50 us time resolution
crop_start_s = float(np.loadtxt(sys.argv[3]))

t_start_MJD = lastMJD + crop_start_s / 86400
dt_MJD = 50e-6 / 86400
t = (np.arange(I.shape[0])+0.5) * dt_MJD + t_start_MJD

np.savetxt("prof.txt", np.vstack([t, I]).T, fmt="%.20f")