#!/usr/bin/env python
import argparse
import getpass
import glob
import math
import os
import sys

# Global constants
NCODIFPARALLEL = 8  # Number of vcraft conversions to do at a time
VCRAFTTIME = "02:00"  # Time to request for vcraft conversion
VCRAFTMEM = 200  # Memory to request per cpu for vcraft conversion


def _main():
    args = get_args()

    keepcodif = args.keep  # Don't rerun CRAFTConverter

    correlateseconds = 20
    if args.bits == 4:
        correlateseconds = 6
    elif args.bits == 8:
        correlateseconds = 4
    elif args.bits == 16:
        correlateseconds = 3

    vcraftfiles = find_vcraft(args.fileglob)

    npol = len(vcraftfiles)
    nant = len(vcraftfiles[0])

    startmjd = args.startmjd
    xfreqs = []

    for file in vcraftfiles[0]:
        hdrfile = file + ".hdr"
        thisvals = parse_vcraft_hdr(hdrfile)

        # as we go over all the headers, search for the earliest start
        # Duration of data in MJD
        durmjd = float(thisvals["NSAMPS_REQUEST"]) / (
            float(thisvals["SAMP_RATE"]) * 86400
        )
        thismjd = float(thisvals["TRIGGER_MJD"]) - durmjd
        if thismjd < startmjd:
            startmjd = thismjd

        # get all of the freqs as we go
        xfreqs += thisvals["FREQS"]

    # Make sure we only have unique freqs and they're ordered
    xfreqs = list(set(xfreqs))
    xfreqs.sort()

    # Double check that the frequencies match
    if npol > 1:
        yfreqs = []
        for file in vcraftfiles[1]:
            hdrfile = file + ".hdr"
            thisvals = parse_vcraft_hdr(hdrfile)
            yfreqs += thisvals["FREQS"]

    yfreqs = list(set(xfreqs))
    yfreqs.sort()

    assert xfreqs == yfreqs, (
        "Frequencies not matching between pols!\n"
        f"X: {xfreqs}\n"
        f"Y: {yfreqs}"
    )

    # These should be constant between header files, so can just grab at end
    beamra = float(thisvals["BEAM_RA"])
    beamdec = float(thisvals["BEAM_DEC"])
    mode = int(thisvals["MODE"])

    # round start MJD down to previous integer second
    startmjd = math.floor(startmjd * 60 * 60 * 24) / (60 * 60 * 24)

    rastr, decstr = posradians2string(
        beamra * math.pi / 180, beamdec * math.pi / 180
    )

    # Overwrite the RA/Dec with command line values if supplied
    if args.ra != None:
        rastr = args.ra
    if args.dec != None:
        decstr = args.dec

    # Write the obs.txt file
    write_obs(rastr, decstr, correlateseconds, startmjd)

    # Write the chandefs file
    write_chandefs(xfreqs, npol)

    if args.ts > 0:
        print("Waiting on CRAFTConverter to finish")
        ret = os.system(f"tsp -S {args.ts}")
        if ret != 0:
            sys.exit(ret)

    # Convert vcraft files
    convertlines = []
    for i in range(npol):
        for f in vcraftfiles[i]:
            runline = get_convert_vcraft_cmd(f, i, keepcodif, args.ts)

            # runline is None if we don't need to re-do this file
            if runline:
                # convert now or do with slurm later
                if args.slurm:
                    convertlines.append(runline)
                else:
                    print(runline)
                    ret = os.system(runline)
                    if ret != 0:
                        sys.exit(ret)

    antnames = [f.split("/")[-1].split("_")[0] for f in vcraftfiles[0]]
    antlist = ",".join(antnames)

    if args.slurm and len(convertlines) > 0:
        write_convert_vcraft_script(convertlines, NCODIFPARALLEL)
        convert_vcraft_slurm(NCODIFPARALLEL)

    # Write a machines file and a run.sh file
    write_run(nant)

    if args.ts > 0:
        print("Waiting on CRAFTConverter to finish")
        ret = os.system("tsp -w")
        if ret != 0:
            sys.exit(ret)

    # Print out the askap2difx command line to run
    # (ultimately, could just run here)
    runline = get_askap2difx_cmd(args, antlist, npol)

    print("\nNow run:")
    print(runline)
    with open("runaskap2difx", "w") as runaskap:
        runaskap.write(runline)
    os.chmod("runaskap2difx", 0o775)


def get_args() -> argparse.Namespace:
    """Parse and validate commmand line arguments

    :return: Command line arguments as a Namespace
    :rtype: :class:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--ra", help="Force RA value")
    parser.add_argument(
        "-d",
        "--dec",
        help="Force Dec value: use no space if declination is negative, i.e., -d-63:20:23.3",
    )
    parser.add_argument(
        "-b", "--bits", type=int, default=1, help="Number of bits"
    )
    parser.add_argument(
        "-k",
        "--keep",
        default=False,
        action="store_true",
        help="Keep exisiting codif files",
    )
    parser.add_argument(
        "-f", "--fpga", help="FPGA and card for delay correction. E.g. c4_f0"
    )
    parser.add_argument(
        "-p", "--polyco", help="Bin config file for pulsar gating"
    )
    parser.add_argument(
        "-i", "--integration", type=float, help="Correlation integration time"
    )
    parser.add_argument(
        "--ts",
        default=0,
        type=int,
        help="Use taskspooler to run CRAFTConverter, with N parallel tasks",
    )
    parser.add_argument(
        "-s",
        "--slurm",
        default=False,
        action="store_true",
        help="Use slurm batch jobs rather than running locally",
    )
    parser.add_argument(
        "-n", "--nchan", type=int, help="Number of spectral channels"
    )
    parser.add_argument(
        "--forceFFT",
        default=False,
        action="store_true",
        help="Force FFT size to equal number of channels (don't increase to 128)",
    )
    parser.add_argument(
        "fileglob", help="glob pattern for vcraft files", nargs="+"
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
    parser.add_argument(
        "--calconly",
        default=False,
        action="store_true",
        help="Stop after creating .calc file",
    )
    parser.add_argument("--dir", help="Directory of local difx files")
    parser.add_argument(
        "--startmjd",
        default=-1,
        type=float,
        help="Force a particular start MJD without searching files for one. Optional.",
    )
    args = parser.parse_args()

    # Check that sensible options were given for the queue destination
    if args.large and not args.gstar:
        parser.error(
            "You can't run large if runnning on skylake (the default, i.e. you didn't use --gstar"
        )
    if args.gstar and args.numskylakenodes > 1:
        parser.error(
            "You can't set the number of skylake nodes if you are running on gstar"
        )

    vcraftglobpattern = args.fileglob
    if len(vcraftglobpattern) > 2:
        parser.error(
            "Can only have at most two fileglobs, corresponding to X and Y pols"
        )

    return args


def posradians2string(rarad: float, decrad: float) -> "tuple[str]":
    """Convert RA and Dec from radians to a pair of strings (hms, dms)

    :param rarad: Right ascension in radians
    :type rarad: float
    :param decrad: Declination in radians
    :type decrad: float
    :return: The RA and Dec as a pair of strings in hms and dms format
        respectively.
    :rtype: tuple[str]
    """
    rah = rarad * 12 / math.pi
    rhh = int(rah)
    rmm = int(60 * (rah - rhh))
    rss = 3600 * rah - (3600 * rhh + 60 * rmm)
    decd = decrad * 180 / math.pi
    decformat = "+%02d:%02d:%010.7f"
    if decd < 0:
        decd = -decd
        decformat = "-" + decformat[1:]
    ddd = int(decd)
    dmm = int(60 * (decd - ddd))
    dss = 3600 * decd - (3600 * ddd + 60 * dmm)
    rastring = "%02d:%02d:%011.8f" % (rhh, rmm, rss)
    decstring = decformat % (ddd, dmm, dss)
    return rastring, decstring


def find_vcraft(fileglobs: "list[str]") -> "list[list[str]]":
    """Find vcraft files based on the provided glob strings.

    :param fileglobs: Strings to glob to find vcraft files.
    :type fileglobs: list[str]
    :return: A list of paths to files as strs for every glob string
        provided
    :rtype: list[str[str]]
    """
    vcraftfiles = []
    for g in fileglobs:
        vcraftfiles.append(glob.glob(g))
        if len(vcraftfiles[-1]) == 0:
            raise FileNotFoundError("Didn't find any vcraft files!")
    if not len(vcraftfiles[0]) == len(vcraftfiles[-1]):
        raise Exception("Number of vcraft files for X and Y doesn't match")

    return vcraftfiles


def parse_vcraft_hdr(hdrfile: str) -> dict:
    """Parse provided vcraft header file and return it as a dictionary.

    :param hdrfile: File to parse
    :type hdrfile: str
    :return: Dictionary of fields in the header file
    :rtype: dict
    """
    vals = {}
    with open(hdrfile) as hdr:
        for line in hdr:
            # Cut off comment from line
            line = line[: line.find("#")]

            # Split line into a list with field name as first value
            line = line.split()

            if len(line) > 2:
                vals[line[0]] = line[1:]
            else:
                vals[line[0]] = line[1]

    # special field: FREQS
    print(vals["FREQS"])
    vals["FREQS"] = vals["FREQS"].split(",")

    return vals


def write_obs(
    rastr: str, decstr: str, correlateseconds: int, startmjd: float
) -> None:
    """Write the obs.txt file containing the start and stop MJD of the
    data and the position to correlate about

    :param rastr: Right ascension as an hh:mm:ss string
    :type rastr: str
    :param decstr: Declination as a dd:mm:ss string
    :type decstr: str
    :param correlateseconds: Duration of the correlation as determined
    by the number of bits per data point
    :type correlateseconds: int
    :param startmjd: Start time of the data as MJD
    :type startmjd: float
    """
    stopmjd = startmjd + correlateseconds / 86400
    output = open("obs.txt", "w")
    output.write(f"startmjd    = {startmjd:.9f}\n")
    output.write(f"stopmjd     = {stopmjd:.9f}\n")
    output.write("srcname     = CRAFTSRC\n")
    output.write(f"srcra       = {rastr}\n")
    output.write(f"srcdec      = {decstr}\n")
    output.close()


def write_chandefs(freqs: "list[str]", npol: int) -> None:
    """Write the chandefs file containing channel definitions. Currently
    vcraft headers have a 1 MHz frequency offset - this is corrected
    for here.

    :param freqs: List of frequencies in MHz as strings as determined
        from the header files
    :type freqs: list[str]
    :param npol: Number of polarisations being processed
    :type npol: int
    """
    output = open("chandefs.txt", "w")
    for i in range(npol):
        for f in freqs:
            print(i, f)
            # Correct the 1 MHz offset in the headers
            # WARN This should probably be regularly checked!
            # Also this can be upper sideband in some cases!
            output.write("%s L 1.185185185185185185\n" % str(int(f) - 1))

    output.close()


def get_convert_vcraft_cmd(
    vcraft: str, polidx: int, keepcodif: bool, ts: int
) -> str:
    """Get the command that will convert a vcraft file into codif

    :param vcraft: vcraft file to be converted
    :type vcraft: str
    :param polidx: Index of the polarisation for this file
    :type polidx: int
    :param keepcodif: If True, only return the command if the codif file
        doesn't already exist (returning None if it does). Otherwise,
        always re-run conversion.
    :type keepcodif: bool
    :param ts: If greater than zero, include tsp call in runline
    :type ts: int
    :return: Command (including arguments) that will convert the vcraft
        file into codif
    :rtype: str
    """
    if not os.path.exists(".bat0"):
        ret = os.system("bat0.pl %s" % (vcraft))
        if ret != 0:
            sys.exit(ret)

    antname = vcraft.split("/")[-1].split("_")[0]
    if antname == "":
        print("Didn't find the antenna name in the header!")
        sys.exit()

    codifname = f"{antname}.p{polidx}.codif"

    # Always run CRAFTConverter if the codif file doesn't exist yet, or
    # if we're forcing a re-run by keepcodif = False
    if not keepcodif or not os.path.exists(codifname):
        runline = f"CRAFTConverter {vcraft} {codifname}"
        if ts > 0:
            runline = "tsp " + runline
    else:
        runline = None

    return runline


def write_convert_vcraft_script(
    convertlines: "list[str]", ncodifparallel: int
) -> None:
    """Write miniscripts that will be used to convert vcraft files to
    codif with slurm

    :param convertlines: Commands (including arguments) that will
        convert vcraft files to codif
    :type convertlines: list[str]
    :param ncodifparallel: Number of conversions to be done in parallel
    :type ncodifparallel: int
    """
    currentuser = getpass.getuser()

    # initialise all the convertcodif miniscripts
    for i in range(ncodifparallel):
        output = open("convertcodif.%d" % (i + 1), "w")
        output.write("#!/bin/bash\n")
        # TODO: change this using home directory!
        output.write(f". /home/{currentuser}/setup_difx\n")
        output.close()
        os.system("chmod 775 convertcodif.%d" % (i + 1))

    # spread out the conversion commands evenly across the miniscripts
    for count, runline in enumerate(convertlines):
        output = open("convertcodif.%d" % ((count % ncodifparallel) + 1), "a")
        output.write(runline + "\n")
        output.close()


def convert_vcraft_slurm(ncodifparallel: int) -> None:
    """Run an overall sbatch script for the CRAFT Conversion stage.
    This will parallelise over ncodifparallel nodes

    :param ncodifparallel: Number of conversions to be run in parallel
    :type ncodifparallel: int
    """
    output = open("runcraftconversionbatch.sh", "w")
    output.write("#!/bin/bash\n")
    output.write("#\n")
    output.write("#SBATCH --job-name=test_craftconverter\n")
    output.write("#SBATCH --output=testt_craftconverter.txt\n")
    output.write("#\n")
    output.write("#SBATCH --ntasks=1\n")
    output.write(f"#SBATCH --time={VCRAFTTIME}\n")
    output.write(f"#SBATCH --mem-per-cpu={VCRAFTMEM}\n")
    output.write(f"#SBATCH --array 1-{ncodifparallel}\n\n")
    output.write("srun ./convertcodif.$SLURM_ARRAY_TASK_ID\n")
    output.close()

    # Run that sbatch script
    os.system("sbatch --wait runcraftconversionbatch.sh")


def write_run(nant: int) -> None:
    """Write a script `run.sh` to be used in the codif to difx
    conversion

    :param nant: Number of antennas in the dataset
    :type nant: int
    """
    output = open("machines", "w")
    for i in range(nant + 2):
        output.write("localhost\n")
    output.close()

    output = open("run.sh", "w")
    output.write("#!/bin/sh\n\n")
    output.write("rm -rf craft.difx\n")
    output.write("rm -rf log*\n")
    output.write("errormon2 6 &\n")
    output.write("export ERRORMONPID=$!\n")
    output.write(
        "mpirun -machinefile machines -np %d mpifxcorr craft.input\n"
        % (nant + 2)
    )
    output.write("kill $ERRORMONPID\n")
    output.write("rm -f craft.difxlog\n")
    output.write("mv log craft.difxlog\n")
    output.close()


def get_askap2difx_cmd(
    args: argparse.Namespace, antlist: "list[str]", npol: int
) -> str:
    """Determine the command (with arguments) to run the codif to difx
    conversion

    :param args: Command line arguments
    :type args: :class:`argparse.Namespace`
    :param antlist: List of antenna names as strings
    :type antlist: list[str]
    :param npol: Number of polarisations being processed
    :type npol: int
    :return: Command (with arguments) to run codif to difx conversion
    :rtype: str
    """
    framesize = 8256 if args.bits == 4 else 8064

    runline = (
        f"{args.dir}/askap2difx.py fcm.txt obs.txt chandefs.txt "
        f"--ants={antlist[:-1]} "
        f"--bits={args.bits} "
        f"--framesize={framesize} "
        f"--npol={npol}"
    )

    if args.slurm:
        runline += " --slurm"
    if args.fpga is not None:
        runline += f" --fpga {args.fpga}"
    if args.polyco is not None:
        runline += f" --polyco={args.polyco}"
    if args.integration is not None:
        runline += f" --integration={args.integration}"
    if args.nchan is not None:
        runline += f" --nchan={args.nchan}"
    if args.calconly:
        runline += " --calconly"
    if args.forceFFT:
        runline += " --forceFFT"
    if args.gstar:
        runline += " --gstar"
    if args.large:
        runline += " --large"
    if args.numskylakenodes > 1:
        runline += " --numskylakenodes=" + str(args.numskylakenodes)
    runline += "\n"
    return runline


if __name__ == "__main__":
    _main()
