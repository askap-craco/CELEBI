# imports
import numpy as np
import glob, subprocess, argparse

def get_args():
    parser = argparse.ArgumentParser(description="Filter antennas based on expected polarisations")
    parser.add_argument("--pols", nargs='+', required=True, help="List of expected polarisations (e.g., X Y or just X)")
    return parser.parse_args()

args = get_args()
expected_pols = args.pols

# glob all files
in_files = glob.glob("*_*_*_f.npy")
out_files = []

# output diagnostic file
diagnostic_out = open("antenna_filtering.txt", "w")

if len(in_files) == 0:
    diagnostic_out.write("No files found to filter.\n")
    diagnostic_out.close()
    exit(0)

# assuming the format output paths follow {prefix}_{antno}_{pol}_f.npy
antnos = list(set([int(file.split('_')[-3]) for file in in_files]))
filtered_antnos = []

for antno in antnos:
    filter_flag = False
    ant_pol_files = []
    
    # Extract the string prefix used for this antenna's files dynamically
    prefix = [f for f in in_files if f"_{antno}_" in f][0].split(f"_{antno}_")[0]

    for pol in expected_pols:
        expected_file = f"{prefix}_{antno}_{pol}_f.npy"
        ant_pol_files.append(expected_file)
        
        # Condition 1: The file entirely failed to generate
        if expected_file not in in_files:
            filter_flag = True
            diagnostic_out.write(f"{expected_file}".ljust(50) + " MISSING\n")
        else:
            # Condition 2: The file generated but has zero data
            farr = np.load(expected_file, mmap_mode = 'r')
            if farr.size == 0:
                filter_flag = True
                diagnostic_out.write(f"{expected_file}".ljust(50) + " ZERO-SIZE\n")
            else:
                diagnostic_out.write(f"{expected_file}".ljust(50) + " OK\n")
        
    diagnostic_out.write("\n")

    # Only pass the antenna if ALL explicitly requested polarisations exist and are non-zero
    if not filter_flag:
        filtered_antnos.append(antno)
        out_files.extend(ant_pol_files)
        
        # create copies of symbolic link as output to nextflow process
        for file in ant_pol_files:
            subprocess.run(f"cp -P {file} {file[:-4]+'_filtered.npy'}", shell=True)

# final diagnostic printing
diagnostic_out.write(f"Unfiltered antennas:  {len(antnos)}\n")
diagnostic_out.write(", ".join([str(a) for a in antnos]) + "\n\n")

diagnostic_out.write(f"Filtered antennas (Kept):  {len(filtered_antnos)}/{len(antnos)}\n")
diagnostic_out.write(", ".join([str(a) for a in filtered_antnos]) + "\n")

diagnostic_out.close()
