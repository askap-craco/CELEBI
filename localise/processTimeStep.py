"""Performs a correlation on a single card-FPGA pair's voltages
"""

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
    if args.polyco:
        polyco = os.path.abspath(args.polyco)
    else:
        polyco = None
    fcm = os.path.abspath(args.fcm)
    topdir = os.getcwd()

    outdir = args.o
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)

    beamdirs = find_vcraft(datadir, beam=args.beam, card=args.card)

    # select specific card/FPGA to process and process it
    freqlabel = args.freqlabel
    print(f"Going to process {freqlabel}")
    if not os.path.exists(freqlabel):
        os.mkdir(freqlabel)
    os.chdir(freqlabel)
    subprocess.call(f"cp {fcm} fcm.txt", shell=True)

    v2oargs = create_v2oargs(args, polyco, datadir, beamdirs)
    v2ocmd = f"python3 {args.dir}/vcraft2obs.py {v2oargs}"

    ret = runCommand(v2ocmd, "vcraft2obs.log")
    if ret != 0:
        print("vcraft2obs failed! (", ret, ")")
        sys.exit(ret)

    find_eop(topdir)

    # Convert CODIF data to difx format
    ret = runCommand("./runaskap2difx", "askap2difx.log")
    if ret != 0:
        print("askap2difx failed! (", ret, ")")
        sys.exit(ret)

    # if .bat0 does not exist in upper directory, copy back
    if not os.path.exists("../../.bat0") and os.path.exists(".bat0"):
        subprocess.call("cp .bat0 ../..", shell=True)

    # If we only need to calculate the interferometer model, quit now
    if args.calconly:
        sys.exit()

    # Launch/process final job
    if args.slurm:
        print("Launching job!")
        subprocess.call("./launchjob", shell=True)
    else:
        subprocess.call("./run.sh", shell=True)
        if args.ref is not None:
            subprocess.call("./run_fill_DiFX", shell=True)
        else:
            subprocess.call("./run_tscrunch_DiFX", shell=True)
        subprocess.call("./runmergedifx", shell=True)
        if args.correctfpgadelays:
            subprocess.call("findOffsets.py", shell=True)
            subprocess.call("./run.sh", shell=True)
            if args.ref is not None:
                subprocess.call("./run_fill_DiFX", shell=True)
            else:
                subprocess.call("./run_tscrunch_DiFX", shell=True)
            subprocess.call("rm -rf craftfrbD2D*", shell=True)
            subprocess.call("./runmergedifx", shell=True)

    print("Done with job")
    os.chdir("../")
    subprocess.call("pkill difxlog", shell=True)
    sys.exit()


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
    parser.add_argument(
        "--uppersideband",
        default=False,
        action="store_true",
        help="Force upper sideband for all channels"
    )
    parser.add_argument(
        "--ref", help="Reference correlation directory", default=None
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

    if args.snoopylog and not os.path.exists(args.snoopylog):
        parser.error(f"Snoopy log file {args.snoopylog} doesn't exist")

    if not os.path.exists(args.fcm):
        parser.error(f"{args.fcm} doesn't exist")

    if args.polyco and not os.path.exists(args.polyco):
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

    # For each antenna subdir, get the beam subdir(s) and look for vcraft
    # files inside of them, stopping once we've verified that any vcraft
    # files exist.

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


def create_v2oargs(
    args: argparse.Namespace,
    polyco: str,
    datadir: str,
    beamdirs: "list[str]",
) -> str:
    """Determine arguments to be passed to vcraft2obs based on command
    line arguments.

    :param args: Command line arguments
    :type args: :class:`argparse.Namespace`
    :param polyco: Absolute path to polyco
    :type polyco: str
    :param datadir: Absolute path to data directory
    :type datadir: str
    :param beamdirs: List of paths to directories for each beam being
        processed
    :type beamdirs: list[str]
    :return: String containing arguments to be passed to vcraft2obs.py
    :rtype: str
    """
    v2oargs = f"--dir={args.dir} --startmjd={args.startmjd}"
    if args.keep:
        v2oargs += " -k"
    if args.ra is not None:
        v2oargs += f" -r{args.ra}"
    if args.dec is not None:
        v2oargs += f" -d{args.dec}"
    if not args.bits == "":
        v2oargs += f" --bits={args.bits}"
    if polyco is not None:
        v2oargs += f" --polyco {polyco}"
    if args.integration is not None:
        v2oargs += f" --integration={args.integration}"
    if args.ts is not None:
        v2oargs += f" --ts={args.ts}"
    if args.nchan is not None:
        v2oargs += f" --nchan={args.nchan}"
    if args.forceFFT:
        v2oargs += " --forceFFT"
    if args.gstar:
        v2oargs += " --gstar"
    if args.large:
        v2oargs += " --large"
    if args.numskylakenodes > 1:
        v2oargs += f" --numskylakenodes={args.numskylakenodes}"
    if args.slurm:
        v2oargs += " --slurm"
        homedir = f"{os.path.expanduser('~')}/"
        if (
            os.path.exists(f"{homedir}.eops.new")
            and os.path.getsize(f"{homedir}.eops.new") > 0
        ):
            subprocess.call(f"mv -f {homedir}.eops.new {homedir}.eops".split())
    if args.calconly:
        v2oargs += " --calconly"
    if args.uppersideband:
        v2oargs += " --uppersideband"
    if args.ref is not None:
        v2oargs += f" --ref={args.ref}"

    v2oargs += f" --fpga {args.freqlabel}"

    for beamdir in beamdirs:
        beamname = os.path.basename(beamdir)
        v2oargs += f' "{datadir}/ak*/{beamname}/*{args.freqlabel}*vcraft"'

    return v2oargs


def runCommand(command: str, log: str) -> int:
    """Run a command, save output in a log file, and return exit code

    :param command: Command to be executed (including arguments)
    :type command: str
    :param log: File to write command STDOUT + STDERR to
    :type log: str
    :return: Exit code of the commmand
    :rtype: int
    """
    print(command)
    proc = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    with open(log, "w") as log_file:
        outs, errs = proc.communicate()
        print(outs.decode("utf-8"))
        print(errs.decode("utf-8"))
        log_file.write(outs.decode("utf-8"))
        log_file.write(errs.decode("utf-8"))

    return proc.returncode


def find_eop(topdir: str) -> None:
    """Find and, if necessary, create Earth Orientation Parameters file

    :param topdir: Directory this script was run from
    :type topdir: str
    """
    if not os.path.exists("eop.txt"):
        topEOP = f"{topdir}/eop.txt"
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
                ret = subprocess.run(f"getEOP.py -l {mjd} > {topEOP}", shell=True).returncode
                if ret != 0:
                    print(
                        "WARNING: getEOP.py call not successful. Your eop.txt "
                        "file is probably empty"
                    )
                    sys.exit(ret)
            else:
                print("Could not find MJD in obs.txt")
                sys.exit()

        print("Copying EOP from top dir")
        print(f"cp {topEOP} eop.txt")
        subprocess.run(f"cp {topEOP} eop.txt", shell=True)


if __name__ == "__main__":
    _main()
