# imports
import numpy as np
import glob, subprocess


# glob all files and find which are size of zero
in_files = glob.glob("*.npy")
out_files = []

# output diagnostic file
diagnostic_out = open("antenna_filtering.txt", "w")

# we have to remove all antenna pols that have size of zero, however,
# in may be the case that only 1 pol of 1 or more antennas were flagged, 
# hence we want to flag the other pol. 

# assuming the format output paths follow *_{antno}_{pol}_f.npy
antnos = [int(file.split('_')[-3]) for file in in_files]
antnos = list(set(antnos))      # make sure not doubling up
filtered_antnos = []

# loop through each antenna and check if both polarisations are there, we can 
# do a shortcut by knowning only the X and Y pols should be here, so just check if
# 2 files exist within out_files
for antno in antnos:
    filter_flag = False
    pol_files = glob.glob(f"*_{antno}_*_f.npy")

    # check size of files
    zero_files = []
    for file in pol_files:
        farr = np.load(file, mmap_mode = 'r')
        if farr.size == 0:
            filter_flag = True
            zero_files += [file]

            # output diagnostics to txt file
            diagnostic_out.write(f"{file}".ljust(50) + " NO\n")
        else:
            diagnostic_out.write(f"{file}".ljust(50) + " YES\n")
        
    diagnostic_out.write("\n")

    if not filter_flag:
        filtered_antnos += [antno]


    if len(zero_files) == 0:
        # both X and Y must exist, hence add them to the filtered antenna pool
        out_files += pol_files

        # create copies of symbolic link as output to nextflow process
        [subprocess.run(f"cp -P {file} {file[:-4]+'_filtered.npy'}", shell = True) for file in pol_files]


# final diagnostic printing
diagnostic_out.write(f"Unfiltered antennas:  {len(antnos)}/{len(antnos)}\n")
[diagnostic_out.write(f"{antno}, ") for antno in antnos]
diagnostic_out.write("\n\n")
diagnostic_out.write(f"Filtered antennas:  {len(filtered_antnos)}/{len(antnos)}\n")
[diagnostic_out.write(f"{filtered_antno}, ") for filtered_antno in filtered_antnos]

