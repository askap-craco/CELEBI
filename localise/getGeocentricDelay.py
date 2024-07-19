import glob
import math
import os
import sys

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

askap_lat = -26.697
askap_lon = 116.631
askap_height = 361
askap_radius = 6374217  # distance from geocentre
c = 299792458.0


def _main():
    rawdata, snoopy_file, numfinderbins, searchms = get_args()
    hdrfiles = find_hdrs(rawdata)

    minfreq, frbpos, triggermjd, samprate, nsamps = parse_hdrs(hdrfiles)

    geocentricdelay = calc_geocentric_delay(triggermjd, frbpos)
    corrstartmjd = calc_corr_start(triggermjd, nsamps, samprate)

    print(f"Geocentric delay is {geocentricdelay}")
    print(f"Lowest frequency is {minfreq}")

    # Determine and print snoopylog2frbgate.py command to be run
    print(get_sl2fg_cmd(minfreq, geocentricdelay, corrstartmjd, snoopy_file, numfinderbins, searchms))


def get_args() -> "tuple[str, str]":
    """Parse command line arguments and check there's the right amount.

    Also verifies that the provided snoopy_file exists, and exits if it
    doesn't.

    :return: rawdata, snoopy_file paths as given in command line
        arguments
    :rtype: tuple[str, str]
    """
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print(f"Usage: {sys.argv[0]} <rawdata directory> <snoopy file> [numfinderbins=7] [searchms=70]")
        print(
            "Raw data directory should contain akXX/beamXX/*.vcraft.hdr files"
        )
        sys.exit()

    rawdata = sys.argv[1]
    snoopy_file = sys.argv[2]

    if not os.path.exists(snoopy_file):
        print(f"{snoopy_file} doesn't exist")
        sys.exit()

    numfinderbins = 7
    if len(sys.argv) > 3:
        numfinderbins = int(sys.argv[3])

    searchms = 70
    if len(sys.argv) > 4:
        searchms = float(sys.argv[4])

    return rawdata, snoopy_file, numfinderbins, searchms


def find_hdrs(rawdata: str) -> "list[str]":
    """Search the rawdata directory for all the data header files (for a
    single polarisation in an antenna).

    While searching, check that we have antenna directories, and then
    at least one polarisation subdirectory within each of those antenna
    directories.

    The information required from the header files (lowest frequency and
    approximate FRB direction) aren't antenna or polarisation-dependent,
    so we only need the headers for a single polarisation of a single
    antenna.

    :param rawdata: Path to the base of the raw data directory (i.e. the
        directory containing the ak?? antenna sub-directories)
    :type rawdata: str
    :return: List of paths to all header files found
    :rtype: list[str]
    """
    antennas = sorted(glob.glob(f"{rawdata}/ak??/"))
    if len(antennas) == 0:
        print("No antennas found!")
        sys.exit()

    beams = glob.glob(f"{antennas[0]}/beam??/")
    if len(beams) == 0:
        print(f"No beams found in {antennas[0]} (in /beam??/)")
        sys.exit()

    hdrfiles = glob.glob(beams[0] + "/*hdr")
    if len(hdrfiles) == 0:
        print("No vcraft header files found!")
        sys.exit()

    return hdrfiles


def parse_hdrs(
    hdrfiles: "list[str]",
) -> "tuple[float, SkyCoord, float, float, int]":
    """Parse vcraft headers to determine the lowest frequency,
    approximate FRB position, MJD of trigger, sample rate, and number of
    samples.

    :param hdrfiles: List of paths to all header files for a single
        polarisation of a single antenna
    :type hdrfiles: list[str]
    :return: minfreq, frbpos, triggermjd, samprate, nsamps
        minfreq: the minimum frequency found (in MHz)
        frbpos: the approximate position of the FRB as determined by the
            central RA and Dec of the beam.
        triggermjd: MJD of the trigger for the voltage dump
        samprate: sample rate per second
        nsamps: number of samples in data
    :rtype: tuple[float, :class:`SkyCoord`, float, float, int]
    """
    minfreq = 99999999
    for hf in hdrfiles:
        with open(hf) as headerin:
            lines = headerin.readlines()
            for line in lines:
                if "FREQ" in line:
                    freqs = line.split()[1].split(",")
                    for f in freqs:
                        if float(f) < minfreq:
                            minfreq = float(f)
                if "BEAM_RA" in line:
                    beamra = float(line.split()[1])
                if "BEAM_DEC" in line:
                    beamdec = float(line.split()[1])
                if "TRIGGER_MJD" in line:
                    triggermjd = float(line.split()[1])
                if "SAMP_RATE" in line:
                    samprate = float(line.split()[1])
                if "NSAMPS_REQUEST" in line:
                    nsamps = int(line.split()[1])

    frbpos = SkyCoord(beamra, beamdec, unit="deg")

    minfreq -= 1    # headers are off by 1 MHz

    return minfreq, frbpos, triggermjd, samprate, nsamps


def calc_geocentric_delay(triggermjd: float, frbpos: SkyCoord) -> float:
    """Calculate the geocentric delay based on the time of the trigger
    and approximate position of the FRB.

    :param triggermjd: MJD of the trigger for the voltage dump
    :type triggermjd: float
    :param frbpos: Approximate position of the FRB (the central RA and
        Dec of the beam)
    :type frbpos: :class:`SkyCoord`
    :return: Geocentric delay in seconds
    :rtype: float
    """
    askap = EarthLocation(
        lon=askap_lon * u.deg, lat=askap_lat * u.deg, height=askap_height * u.m
    )
    t = Time(triggermjd, format="mjd")
    altaz = frbpos.transform_to(AltAz(obstime=t, location=askap))
    geocentricdelay = math.sin(altaz.alt.radian) * askap_radius / c

    """
    Would be much simpler to just get the elevation from the vcraft 
    file!! But this is for the centre of the PAF, so could be a bit off 
    vs the elevation of the actual beam (which is itself already a bit 
    off vs the actual FRB, but within half a degree which is good 
    enough)
    geocentricdelay = math.sin(ant_el_rad) * askap_radius / c
    """

    return geocentricdelay


def calc_corr_start(triggermjd: float, nsamps: int, samprate: float) -> float:
    """Calculate the start time to use for the correlation in MJD.

    This is calculated as the latest whole second before the start of
    the dumped data, i.e.
        `t_start = t_trigger - n_samps/samprate`
    rounded down to the nearest whole second.

    :param triggermjd: MJD of the time of the trigger
    :type triggermjd: float
    :param nsamps: Number of samples in the data
    :type nsamps: int
    :param samprate: Sample rate of the data in samples per second
    :type samprate: float
    :return: MJD of the calculated correlation start time
    :rtype: float
    """
    corrstartmjd = math.floor(triggermjd * 86400 - nsamps / samprate) / 86400.0
    return corrstartmjd


def get_sl2fg_cmd(
    minfreq: float,
    geocentricdelay: float,
    corrstartmjd: float,
    snoopy_file: str,
    numfinderbins: int,
    searchms: float,
) -> str:
    """Determine the snoopylog2frbgate.py command to be run next.

    :param minfreq: Minimum frequency of data (in MHz)
    :type minfreq: float
    :param geocentricdelay: Geocentric delay (in s)
    :type geocentricdelay: float
    :param corrstartmjd: Correlation start time (in MJD)
    :type corrstartmjd: float
    :param snoopy_file: Snoopy candidate file
    :type snoopy_file: str
    :param numfinderbins: Number of finder bins to run
    :type numfinderbins: int
    :param searchms: Search space to run finder bins over
    :type searchms: float
    :return: snoopylog2frbgate.py command call
    :rtype: str
    """
    sl2fgcmd = "snoopylog2frbgate.py"
    sl2fgcmd += f" -f {minfreq:.3f}"
    sl2fgcmd += f" --timediff {geocentricdelay * 1e3:.3f}"
    sl2fgcmd += f" --corrstartmjd {corrstartmjd:.9f}"
    sl2fgcmd += f" --numfinderbins {numfinderbins}"
    sl2fgcmd += f" --searchms {searchms:.3f}"
    sl2fgcmd += f" {snoopy_file}"
    return sl2fgcmd


if __name__ == "__main__":
    _main()
