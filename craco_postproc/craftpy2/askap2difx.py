#!/usr/bin/env python
import argparse
import getpass
import math
import os
import sys

## Argument parser
parser = argparse.ArgumentParser()
parser.add_argument("fcm", help="ASKAP .fcm file describing array")
parser.add_argument(
    "obs",
    help="Flat text .obs file containing start/stop time, source position, baseband files",
)
parser.add_argument(
    "chan",
    help="Flat text file containing 1 line per subband, centre freq, sideband, and bandwidth",
)
parser.add_argument(
    "-a",
    "--ants",
    help="Comma separated list of antenna names e.g., ak01,ak02,ak03 etc.  All must be present in .fcm, and all must have a akxx.codif file present in this directory.",
)
parser.add_argument(
    "-b", "--bits", help="Number of bits/sample. Default 1", type=int
)
parser.add_argument(
    "-t", "--threads", help="Number of DIFX threads. Default 8", type=int
)
parser.add_argument(
    "-F", "--framesize", help="Codif framesize. Default 8064", type=int
)
parser.add_argument(
    "-f", "--fpga", help="FPGA and card for delay correction. E.g. c4_f0"
)
parser.add_argument("-p", "--polyco", help="Bin config file for pulsar gating")
parser.add_argument(
    "-i",
    "--integration",
    default="1.3824",
    help="Correlation integration time",
)
parser.add_argument(
    "-n", "--nchan", type=int, default=128, help="Number of spectral channels"
)
parser.add_argument(
    "-s",
    "--slurm",
    default=False,
    action="store_true",
    help="Use slurm batch jobs rather than running locally",
)
parser.add_argument(
    "--forceFFT",
    default=False,
    action="store_true",
    help="Force FFT size to equal number of channels (don't increase to 128)",
)
parser.add_argument(
    "--npol",
    help="Number of polarisations",
    type=int,
    choices=[1, 2],
    default=1,
)
parser.add_argument(
    "--gstar",
    default=False,
    action="store_true",
    help="Set if using gstar for correlation",
)
parser.add_argument(
    "--large",
    default=False,
    action="store_true",
    help="Set if 32 nodes, 384 tasks are required (i.e., 23GB memory needed per task; else 16 nodes, 192 tasks will be used for 11.5GB per task",
)
parser.add_argument(
    "--numskylakenodes", default=1, type=int, help="Use 32x this many CPUs"
)
args = parser.parse_args()


## Replace Mark5B with VDIF in vex file
ret = os.system("sed -i 's/MARK5B/VDIF/g' craftfrb.vex")
if ret != 0:
    exit(1)


## Run updateFreqs
runline = f"updatefreqs.py craftfrb.vex --npol={args.npol} {args.chan}"
if args.nchan is not None:
    runline += f" --nchan={args.nchan}"
print("Running: " + runline)
ret = os.system(runline)
if ret != 0:
    sys.exit(1)

## Update the vex file to say "CODIF" rather than "VDIF"
ret = os.system(f"sed -i -e 's/VDIF5032/CODIFD{framesize}/g' craftfrb.vex")
if ret != 0:
    sys.exit(1)


## Create the machines and threads file and a corresponding run script, or the slurm batch job, as appropriate
if args.slurm:
    currentdir = os.getcwd()
    numprocesses = args.npol * len(datafilelist[0]) + 5
    print(("Approx number of processes", numprocesses))
    numprocessingnodes = numprocesses - (args.npol * len(datafilelist[0]) + 1)

    # Create a run script that will actually be executed by sbatch
    batchout = open("runparallel", "w")
    batchout.write("#!/bin/bash\n")
    batchout.write("#\n")
    batchout.write(f"#SBATCH --job-name=difx_{basename}\n")
    batchout.write(f"#SBATCH --output={basename}.mpilog\n")
    batchout.write("#\n")
    # Request 32xnumskylakenodes CPUs, so it will fit on numskylakenodes
    batchout.write(f"#SBATCH --ntasks={32*args.numskylakenodes:d}\n")
    batchout.write("#SBATCH --time=8:00\n")
    batchout.write("#SBATCH --cpus-per-task=1\n")
    batchout.write("#SBATCH --mem-per-cpu=4000\n\n")
    batchout.write(
        "difxlog {0} {1}/{0}.difxlog 4 &\n\n".format(basename, currentdir)
    )
    batchout.write(
        f"srun -n{numprocesses} --overcommit mpifxcorr {basename}.input --nocommandthread\n"
    )
    batchout.write("./runmergedifx\n")
    batchout.close()

    # Write the threads file
    threadsout = open(f"{basename}.threads", "w")
    threadsout.write(f"NUMBER OF CORES:    {numprocessingnodes}\n")
    for i in range(numprocessingnodes):
        threadsout.write(f"{difxThreads}\n")
    threadsout.close()
else:
    # We just want Ndatastream processes, plus a head node, plus one computational node
    numprocesses = args.npol * len(datafilelist[0]) + 2

    # Write the machines file (all running on localhost)
    machinesout = open("machines", "w")
    for i in range(numprocesses):
        machinesout.write("localhost\n")
    machinesout.close()

    # Write the threads file
    threadsout = open("craftfrb.threads", "w")
    threadsout.write("NUMBER OF CORES:    1\n")
    threadsout.write(f"{difxThreads}\n")
    threadsout.close()

    ## Create a little run file for running the observations
    ## This is used in preference to startdifx because startdifx uses some MPI options we don't want
    runout = open("run.sh", "w")
    runout.write("#!/bin/sh\n\n")
    runout.write("rm -rf craftfrb.difx\n")
    runout.write("rm -rf log*\n")
    runout.write("errormon2 6 &\n")
    runout.write("export ERRORMONPID=$!\n")
    runout.write(
        "mpirun -machinefile machines -np %d mpifxcorr craftfrb.input\n"
        % (numprocesses)
    )
    runout.write("kill $ERRORMONPID\n")
    runout.write("rm -f craftfrb.difxlog\n")
    runout.write("mv log craftfrb.difxlog\n")
    runout.close()
    os.chmod("run.sh", 0o775)

    print("# First run the correlation:")
    runline = "./run.sh\n"
    print(runline)

## Print the line needed to run the stitching and then subsequently difx2fits
print("Then run difx2fits")
runline = "rm -rf {0}D2D.* ; mergeOversampledDiFX.py craftfrb.stitchconfig {0}.difx\n".format(
    basename
)
print(runline)
with open("runmergedifx", "w") as runmerge:
    runmerge.write(runline)
os.chmod("runmergedifx", 0o775)

## Produce the little script needed to run difx2fits for the solveFPGA script
runline = f"difx2fits {basename}D2D\n"
with open("runfpgadifx2fits", "w") as rund2f:
    rund2f.write(runline)
os.chmod("runfpgadifx2fits", 0o775)
