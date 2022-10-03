#!/usr/bin/env python
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# Subroutine to zero noise in the profile
def zeronoise(profile):
    numprofilebins = profile.shape[0]
    profilemin = np.min(profile)
    noiseonly = []
    for val in profile:
        if val < -profilemin:
            noiseonly.append(val)
    stddev = np.std(noiseonly)

    #Zero all the noisy bits of the image
    windowsize = 8
    for i in range(numprofilebins):
        if not np.mean(profile[i-windowsize//2:i+windowsize//2]) > \
               4*stddev/np.sqrt(windowsize):
            profile[i] = 0.0

    #Go back and check for isolated noisy bits, squash them too
    for i in range(numprofilebins):
        if profile[i-1] == 0 and profile[(i+1)%numprofilebins] == 0:
            if profile[i-2] == 0 or profile[(i+2)%numprofilebins] == 0:
                profile[i] = 0.0

#Subroutine to calculate number of bins
def getNumMergedBins(profile, numprofilebins):
    numbins = 0
    for i in range(numprofilebins):
        if profile[i] != profile[i-1]:
            numbins = numbins + 1
    return numbins

#Subroutine to merge smallest delta(profile)
def mergeSmallestDifference(profile, numprofilebins):
    smallestchange = 9e99
    index = 0
    for i in range(numprofilebins):
        deltaprofile = abs(profile[i] - profile[i-1])
        if deltaprofile > 0:
            leftwidth = 1
            rightwidth = 1
            while profile[i-leftwidth-1] == profile[i-1]:
                leftwidth = leftwidth + 1
            while profile[(i+rightwidth)%numprofilebins] == profile[i]:
                rightwidth = rightwidth + 1
            if deltaprofile < smallestchange:
                smallestchange = deltaprofile
                savedlw = leftwidth
                savedrw = rightwidth
                index = i

    if profile[index] == 0.0 or profile[index-1] == 0.0:
        profile[index] = 0.0
        profile[index-1] = 0.0
        return
    
    meanvalue = (savedlw*profile[index-1] + savedrw*profile[index])\
                /(savedlw+savedrw)
    for i in range(savedlw+savedrw):
        profile[(index-savedlw+i)%numprofilebins] = meanvalue
    return

# Main code
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: {0} <timeseries file> <polyco file>".format(sys.argv[0]))
        sys.exit()

    timeseriesfile = sys.argv[1]
    polycofile = sys.argv[2]
    numbins = 20

    if not os.path.exists(timeseriesfile):
        print("Time series file", timeseriesfile, "doesn't exist")
        sys.exit()
    try:
        timeseries = np.load(sys.argv[1])
    except ValueError:
        timeseries = np.transpose(np.loadtxt(sys.argv[1]))
    print(timeseries.shape)
    numprofilebins = timeseries.shape[1]
    binwidthmjd = timeseries[0][1] - timeseries[0][0]

    if not os.path.exists(polycofile):
        print("Polyco file", polycofile, "doesn't exist")
        sys.exit()
    polycolines = open(polycofile).readlines()
    polycomjd = float(polycolines[0].split()[3]) 
    polycoperiod = 1.0/float(polycolines[1].split()[1])

    # Zero the noise in the profile
    zeroedprofile = np.copy(timeseries[1])
    zeronoise(zeroedprofile)

    #Now merge bins until you reach required number
    outputbins = getNumMergedBins(zeroedprofile, numprofilebins)
    print("Number of merged bins is " + str(outputbins))
    while outputbins > numbins:
        mergeSmallestDifference(zeroedprofile, numprofilebins)
        outputbins = getNumMergedBins(zeroedprofile, numprofilebins)

    #Work out the bin phase edges and the values
    binphases = []
    binvalues = []
    lastvalue = zeroedprofile[-1]
    for i in range(numprofilebins):
        if zeroedprofile[i] != lastvalue:
            binphases.append((timeseries[0][i] + binwidthmjd/2.0 - polycomjd) * 86400 / polycoperiod)
            binvalues.append(lastvalue)
        lastvalue = zeroedprofile[i]
    if len(binphases) != numbins:
        print("Whoops - I somehow averaged down to " + str(len(binphases)) + " bins, rather than " + str(numbins) + " like you asked! Guess we'll use that instead...")
        numbins = len(binphases)

    #Normalise the weights
    weightsum = 0.0
    for i in range(numbins):
        binw = binphases[i] - binphases[i-1]
        if binw < 0.0:
            binw = binw + 1.0
        weightsum = weightsum + binw*binvalues[i]
    for i in range(numbins):
        binvalues[i] = binvalues[i]/weightsum

    #Write the binconfig file
    binconffile = "craftfrb.htrgate.binconfig"
    binconfout  = open(binconffile, 'w')
    binconfout.write("NUM POLYCO FILES:   1\n")
    binconfout.write("POLYCO FILE 0:     {0}\n".format(polycofile))
    binconfout.write("NUM PULSAR BINS:    {0}\n".format(numbins))
    binconfout.write("SCRUNCH OUTPUT:     TRUE\n")
    for i in range(numbins):
        phase = binphases[i]
        binconfout.write(("BIN PHASE END %i:" % (i)).ljust(20))
        binconfout.write("%f\n" % (phase))
        binconfout.write(("BIN WEIGHT %i:"%(i)).ljust(20))
        binconfout.write("%f\n" % (binvalues[i]))

    # Plot both the original profile and the matched filter
    plt.plot(timeseries[0], timeseries[1], 'b-')
    plt.plot(timeseries[0], zeroedprofile, 'r-') #Plot the bin profile
    plt.savefig("craftfrb.matchedfilter.png") 
