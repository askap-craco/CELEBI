import sys
import numpy as np

SNs = sys.argv[1]
inds = sys.argv[2]
#print(list(SNs))
#print(SNs.split())
#print(list(inds))
#print(inds.split())
#print(np.max(float(SNs.split())))

def findMax(SNs, inds):
    #np.argmax(list(SNs))
    #return np.array(list(inds))[np.argmax(list(SNs))]
    inds_split = inds.split()
    return inds_split[np.argmax([float(i) for i in SNs.split()])]
    #return np.array(inds.split())[np.argmax(float(SNs.split()))]

try:
    print(findMax(SNs,inds))
except:
    exit()
    

