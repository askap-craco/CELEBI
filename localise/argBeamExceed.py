import sys
import numpy as np

BMINs = np.array([float(i) for i in sys.argv[1].split()])
BMAXs = np.array([float(i) for i in sys.argv[2].split()])
beamBMIN = float(sys.argv[3])
beamBMAX = float(sys.argv[4])
inds = sys.argv[5]

#print(BMINs, BMAXs, beamBMIN, beamBMAX, inds)

def findIndi(BMINs, BMAXs, beamBMIN, beamBMAX, inds):
    inds_split = np.array(inds.split())
    return inds_split[(BMINs > 3*beamBMIN) | (BMAXs > 3*beamBMAX)]

try:
    print(*findIndi(BMINs, BMAXs, beamBMIN, beamBMAX, inds))
except:
    exit()
