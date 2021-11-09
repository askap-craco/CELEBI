#!/usr/bin/env python3
import argparse
import glob
import os
import re
import subprocess
import sys


def _main():
    # Parse and verify command line arguments
    args = get_args()

    # A few useful paths
    datadir = os.path.abspath(args.timestep)
    polyco = os.path.abspath(args.polyco)
    fcm = os.path.abspath(args.fcm)
    topDir = os.getcwd()

    bins = get_nbins(polyco)

    outdir = args.o

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)

    beamdirs = find_vcraft(datadir, beam=args.beam)
    npol = len(beamdirs)

    # select specific card/FPGA to process and process it
    freqlabels = []
    freqlabel = args.freqlabel
    freqlabels.append(freqlabel)
    print(f"Going to process {freqlabel}")
    if not os.path.exists(freqlabel):
        os.mkdir(freqlabel)
    os.chdir(freqlabel)
    os.system("cp %s fcm.txt" % fcm)
    if os.path.exists("../../eopjunk.txt"):
        os.system("cp ../../eopjunk.txt .")
    # Copy .bat0 file if it exists
    if os.path.exists("../../.bat0"):
        os.system("cp ../../.bat0 .")

    torun = (
        f"{args.dir}/vcraft2obs.py --dir={args.dir} "
        "--startmjd={str(args.startmjd)}"
    )
    if args.keep:
        torun += " -k"
    if args.ra is not None:
        torun = torun + f" -r{args.ra}"
    if args.dec is not None:
        torun = torun + f" -d{args.dec}"
    if not args.bits == "":
        torun = torun + f" --bits={args.bits}"
    if polyco is not None:
        torun += f" --polyco {polyco}"
    if args.integration is not None:
        torun += f" --integration={args.integration}"
    if args.ts is not None:
        torun += f" --ts={args.ts}"
    if args.nchan is not None:
        torun += f" --nchan={args.nchan}"
    if args.forceFFT:
        torun += " --forceFFT"
    if args.gstar:
        torun += " --gstar"
    if args.large:
        torun += " --large"
    if args.numskylakenodes > 1:
        torun += f" --numskylakenodes={args.numskylakenodes}"
    if args.slurm:
        torun += " --slurm"
        homedir = os.path.expanduser("~") + "/"
        if (
            os.path.exists(homedir + ".eops.new")
            and os.path.getsize(homedir + ".eops.new") > 0
        ):
            os.system("mv -f " + homedir + ".eops.new " + homedir + ".eops")
    if args.calconly:
        torun += " --calconly"
    beamname = os.path.basename(beamdirs[0])
    torun += f" --fpga {freqlabel}"
    torun += f'"{datadir}/ak*/{beamname}/*{freqlabel}*vcraft"'
    if npol == 2:
        beamname = os.path.basename(beamdirs[1])
        torun += f' "{datadir}/ak*/{beamname}/*{freqlabel}*vcraft"'

    print(torun)
    ret = runCommand(torun, "vcraft2obs.log")
    if ret != 0:
        print("vcraft2obs failed! (", ret, ")")
        sys.exit(ret)

    # Find and, if necessary, create EOP file
    if not os.path.exists("eop.txt"):
        topEOP = f"{topDir}/eop.txt"
        if not os.path.exists(topEOP):
            mjd = None
            with open("obs.txt") as f:
                for line in f:
                    match = re.search(r"startmjd\s*=\s*(\S+)", line)
                    if match:
                        mjd = match.group(1)
                        break
            if mjd is not None:
                print(f"getEOP.py -l {mjd} > {topEOP}")
                ret = os.system(f"getEOP.py -l {mjd} > {topEOP}")
                if ret != 0:
                    print(
                        "WARNING: getEOP.py call not successful. Your eop.txt file is probably empty"
                    )
                    sys.exit(ret)
            else:
                print("Could not find MJD in obs.txt")
                sys.exit()

        print("Copying EOP from top dir")
        print(f"cp {topEOP} eop.txt")
        os.system(f"cp {topEOP} eop.txt")

    # Convert CODIF data to difx format
    ret = runCommand("./runaskap2difx", "askap2difx.log")
    if ret != 0:
        print("askap2difx failed! (", ret, ")")
        sys.exit(ret)
    # if .bat0 does not exist in upper directory, copy back
    if not os.path.exists("../../.bat0") and os.path.exists(".bat0"):
        os.system("cp .bat0 ../..")

    # If we only need to calculate the interferometer model, quit now
    if args.calconly:
        sys.exit()

    if args.slurm:
        os.system("./launchjob")
    else:
        os.system("./run.sh")
        os.system("./runmergedifx")
        if args.correctfpgadelays:
            os.system("findOffsets.py")
            os.system("./run.sh")
            os.system("rm -rf craftfrbD2D*")
            os.system("./runmergedifx")
    os.chdir("../")


def get_args() -> argparse.Namespace:
    """Parse command line arguments

    :return: Command line argument paramters
    :rtype: :class:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t",
        "--timestep",
        help="Timestep (the directory name to process)",
        required=True,
    )
    parser.add_argument(
        "--name",
        default="CRAFT",
        help="Base name for the output fits files",
    )
    parser.add_argument("-r", "--ra", help="Force RA value")
    parser.add_argument(
        "-d",
        "--dec",
        help="Force Dec value: use no space if declination is "
        "negative, i.e., -d-63:20:23.3",
    )
    parser.add_argument(
        "-b",
        "--bits",
        type=int,
        default=1,
        help="Number of bits. Default 1",
    )
    parser.add_argument(
        "-i",
        "--integration",
        type=float,
        help="Correlation integration time",
    )
    parser.add_argument(
        "-n",
        "--nchan",
        type=int,
        help="Number of spectral channels",
    )
    parser.add_argument(
        "--forceFFT",
        default=False,
        action="store_true",
        help="Force FFT size to equal number of channels (don't increase to "
        "128)",
    )
    parser.add_argument(
        "-f",
        "--fcm",
        default="fcm.txt",
        help="Name of the fcm file",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--polyco",
        help="Bin config file for pulsar gating",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--correctfpgadelays",
        default=False,
        action="store_true",
        help="Figure out and correct 7 microsec FPGA delays",
    )
    parser.add_argument(
        "-S",
        "--suppress",
        default=False,
        action="store_true",
        help="Don't create FITS file",
    )
    parser.add_argument(
        "-B", "--beam", help="Correlate a specific beam: blank means both"
    )
    parser.add_argument(
        "--card",
        default="",
        help="Correlate only a specific card; blank means all",
    )
    parser.add_argument(
        "-k",
        "--keep",
        default=False,
        action="store_true",
        help="Keep existing codif files",
    )
    parser.add_argument(
        "-s",
        "--snoopylog",
        help="Snoopy log file, default blank, if not default "
        "will use this to correlate on-pulse",
        required=True,
    )
    parser.add_argument(
        "--slurm",
        default=False,
        action="store_true",
        help="Use slurm batch jobs rather than running locally",
    )
    parser.add_argument(
        "--ts",
        default=0,
        type=int,
        help="Use taskspooler to run CRAFTConverter, with N " "parallel tasks",
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
        help="Set if 32 nodes, 384 tasks are required (i.e., "
        "23GB memory needed per task; else 16 nodes, 192 "
        "tasks will be used for 11.5GB per task",
    )
    parser.add_argument(
        "--numskylakenodes",
        default=1,
        type=int,
        help="Use 32x this many CPUs",
    )
    parser.add_argument("-o", help="Output directory for data")
    parser.add_argument("--freqlabel", help="Freqlabel to process", type=str)
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
        help="Force a particular start MJD without searching "
        "files for one.",
    )
    args = parser.parse_args()
    verify_args(args, parser)
    return args


def verify_args(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
) -> None:
    """Verify that provided command line arguments are valid

    :param args: Command line arguments
    :type args: :class:`argparse.Namespace`
    :param parser: Parser object used to parse arguments
    :type parser: :class:`argparse.ArgumentParser`
    """

    # Check that sensible options were given for the queue destination
    if args.large and not args.gstar:
        parser.error(
            "You can't run large if runnning on skylake (the default, i.e. you"
            " didn't use --gstar"
        )
    if args.gstar and args.numskylakenodes > 1:
        parser.error(
            "You can't set the number of skylake nodes if you are running on "
            "gstar"
        )

    # Check that provided files exist
    if not os.path.exists(args.timestep):
        parser.error(
            f"Target directory (timestep) {args.timestep} doesn't exist"
        )

    if not os.path.exists(args.snoopylog):
        parser.error(f"Snoopy log file {args.snoopylog} doesn't exist")

    if not os.path.exists(args.fcm):
        parser.error(f"{args.fcm} doesn't exist")

    if not os.path.exists(args.polyco):
        parser.error(f"binconfig file {args.polyco} does not exist")


def get_nbins(polyco: str) -> int:
    """Parse provided bin config file for the number of bins

    :param polyco: Path to bin config file
    :type polyco: str
    :return: Number of bins
    :rtype: int
    """
    # Parse bin config file for number of bins
    nbins = 1
    lines = open(polyco).readlines()
    for i, line in enumerate(lines):
        if "NUM PULSAR BINS" in line:
            nbins = int(line.split(":")[-1].strip())
            if "TRUE" in lines[i + 1]:
                nbins = 1
            break

    return nbins


def find_vcraft(datadir: str, beam: str = None, card: str = "") -> "list[str]":
    """Determine beam directories containing data to be processed.

    Searches within the provided datadir for vcraft files. If none are
    found, exit. If they are found, return a list of paths to the beam
    subdirectories for the first antenna containing data.

    :param datadir: Directory containing the per-antenna directories
    :type datadir: str
    :param beam: Beam subdirectory to return. If None, both beams are
        returned. [Default = None]
    :type beam: str, optional
    :param card: If given, return files for only a specific card.
        [Default = all]
    :type card: str, optional
    :return: List of paths to beam directories to process
    :rtype: list[str]
    """
    examplefiles = []

    # find all antenna subdirectories
    antennadirs = sorted(glob.glob(f"{datadir}/ak*"))

    """
    For each antenna subdir, get the beam subdir(s) and look for vcraft
    files inside of them, stopping once we've verified that any vcraft
    files exist.
    """
    for a in antennadirs:
        if beam is None:
            beamdirs = sorted(glob.glob(f"{a}/*"))
        else:
            beamdirs = [f"{a}/{beam}"]
            if not os.path.exists(f"{a}/{beam}"):
                print(f"{a}/{beam} doesn't exist, aborting")
                sys.exit()

        for b in beamdirs:
            vcraftfiles = glob.glob(f"{b}/*[ac]{card}*vcraft")

            if len(vcraftfiles) > 0:
                examplefiles = sorted(vcraftfiles)
                break
        if len(examplefiles) > 0:
            break

    if len(examplefiles) == 0:
        print("Couldn't find any vcraft files")
        sys.exit()

    npol = len(beamdirs)
    if npol == 0:
        print("Could not find any beams. Aborting")
        sys.exit()

    if npol > 2:
        print("Too many beams found! Aborting")
        sys.exit()

    return beamdirs


def runCommand(command: str, log: str) -> int:
    """Run a command, save output in a log file, and return exit code

    :param command: Command to be executed (including arguments)
    :type command: str
    :param log: File to write command STDOUT + STDERR to
    :type log: str
    :return: Exit code of the commmand
    :rtype: int
    """
    proc = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    with open(log, "w") as log_file:
        outs, errs = proc.communicate()
        print(outs)
        print(errs)
        log_file.write(outs)
        log_file.write(errs)

    return proc.returncode


if __name__ == "__main__":
    _main()
