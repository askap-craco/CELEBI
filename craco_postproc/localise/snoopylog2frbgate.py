#!/usr/bin/env python3
import argparse
import os
import sys

fakepulsarperiod = 10  # seconds


def _main():
    args = get_args()
    validate_args(args)

    timediffsec = args.timediff / 1000.0

    cand = parse_snoopy(args.snoopylog)
    pulsewidthms = float(cand[3]) * 1.7  # filterbank width = 1.7ms
    dm = float(cand[5])
    mjd = float(cand[7])

    # Figure out the time at the midpoint of the pulse
    midmjd = (
        mjd
        - (
            dm * 0.00415 / ((args.freq / 1e3) ** 2)
            - dm * 0.00415 / ((args.freq + 168) / 1e3) ** 2
        )
        / 86400.0
    )

    # Figure out the best int time
    bestinttime = calc_best_int_time(args.corrstartmjd, midmjd)
    print(bestinttime)

    polycorefmjd, hh, mm, ss = calc_polyco_ref_mjd(args.corrstartmjd)

    # Write out the polyco file
    polycopath = write_polyco(polycorefmjd, hh, mm, ss, dm, args.freq)

    # Gate binconfig
    gatebinedges, gateweights = calc_gate_bins(
        cand, timediffsec, polycorefmjd
    )
    write_binconfig(
        "craftfrb.gate.binconfig", polycopath, gatebinedges, gateweights
    )

    # RFI binconfig
    rfibinedges, rfiweights = calc_rfi_bins(gatebinedges)
    write_binconfig(
        "craftfrb.rfi.binconfig", polycopath, rfibinedges, rfiweights
    )

    # High time resolution binconfig
    htrbinedges, htrweights, numbins = calc_htr_bins(cand, gatebinedges[0])
    write_binconfig(
        "craftfrb.bin.binconfig",
        polycopath,
        htrbinedges,
        htrweights,
        scrunch=False,
    )

    # Finder binconfig
    finderbinedges, finderweights, numfinderbins = calc_finder_bins(
        rfibinedges
    )
    write_binconfig(
        "craftfrb.finder.binconfig",
        polycopath,
        finderbinedges,
        finderweights,
        scrunch=False,
    )

    # And write out a little script ready to do the various subtractions
    write_subtractions_script(
        gatebinedges,
        rfibinedges,
        finderbinedges,
        htrbinedges,
        numbins,
        numfinderbins,
    )


def get_args() -> argparse.Namespace:
    """Parse command line arguments

    :return: Command line argument paramters
    :rtype: :class:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser(
        description="Turn a snoopy log into a binconfig and polyco for DiFX."
    )
    parser.add_argument("snoopylog", metavar="S", help="The snoopy log file")
    parser.add_argument(
        "-f",
        "--freq",
        type=float,
        default=-1,
        help="The reference freq at which snoopy DM was calculated",
    )
    parser.add_argument(
        "--timediff",
        type=float,
        default=-99999,
        help="The time difference between the VCRAFT and snoopy log arrival "
        "times for the pulse, including geometric delay, in ms",
    )
    parser.add_argument(
        "--corrstartmjd",
        type=float,
        default=-1,
        help="When the correlation will start",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Validate arguments

    :param args: Command line arguments
    :type args: :class:`argparse.Namespace`
    """
    if args.freq < 0:
        print(
            "You have to supply a frequency for the snoopy files! Don't be "
            "lazy, that's how accidents happen."
        )
        sys.exit()

    if args.timediff < -10000:
        print(
            "You have to specify a timediff! It should just be the geometric "
            "delay from ASKAP to the geocentre, "
        )
        print(
            "now everything has been fixed.  Don't be lazy, that's how "
            "accidents happen."
        )
        sys.exit()

    if args.corrstartmjd < 0:
        print(
            "You have to specify a corrstartmjd. getGeometricDelay.py will "
            "give it to you"
        )
        sys.exit()


def parse_snoopy(snoopy_file: str) -> "list[str]":
    """Parse snoopy file, returning a candidate as a list of strings.

    We expect this file to only contain one candidate (the triggering
    candidate), so we only return one line.

    :param snoopy_file: Path to snoopy candidate file
    :type snoopy_file: str
    :return: Candidate information as a list of strings. Each string is
        a whitespace-separated value in the candidate file.
    :rtype: list[str]
    """
    nocommentlines = []
    for line in open(snoopy_file):
        print(line)
        if len(line) > 1 and not line[0] == "#":
            nocommentlines.append(line)
            print(f"Snoopy info {nocommentlines}")
    if len(nocommentlines) != 1:
        print("ERROR: No information found")
        sys.exit()

    return nocommentlines[0].split()


def calc_best_int_time(corrstartmjd: float, midmjd: float) -> float:
    """Calculate the best integration time such that the FRB should be
    within a single integration

    :param corrstartmjd: Start time of correlation (in MJD)
    :type corrstartmjd: float
    :param midmjd: Central time of the FRB (in MJD)
    :type midmjd: float
    :return: Best integration time (in seconds)
    :rtype: float
    """
    subintsec = 0.13824
    bestinttime = 2 * (midmjd - corrstartmjd) * 86400
    nsubints = int(round(bestinttime / subintsec))
    bestinttime = nsubints * subintsec
    return bestinttime


def calc_polyco_ref_mjd(corrstartmjd: float) -> "tuple[float, int, int, int]":
    """Calculate the start time for the polyco by going back to the
    integer second boundary immediately before the correlation start
    time.

    :param corrstartmjd: Start time of correlation (in MJD)
    :type corrstartmjd: float
    :return: Polyco reference time (in MJD) and the hour, minute, and
        seconds of that time as integers
    :rtype: tuple[float, int, int, int]
    """
    polycorefmjdint = int(corrstartmjd)
    polycorefseconds = int((corrstartmjd - int(corrstartmjd)) * 86400)

    hh = polycorefseconds // 3600
    mm = (polycorefseconds - hh * 3600) // 60
    ss = polycorefseconds - (hh * 3600 + mm * 60)

    polycorefmjd = polycorefmjdint + float(polycorefseconds) / 86400.0

    return polycorefmjd, hh, mm, ss


def write_polyco(
    polycorefmjd: float, hh: int, mm: int, ss: int, dm: float, freq: float
) -> str:
    """Write out the polyco file

    :param polycorefmjd: Polyco reference time (in MJD)
    :type polycorefmjd: float
    :param hh: Hour of the reference time
    :type hh: int
    :param mm: Minute of the reference time
    :type mm: int
    :param ss: Second of the reference time
    :type ss: int
    :param dm: Dispersion measure of the triggering FRB candidate
    :type dm: float
    :param freq: Reference frequency at which snoopy DM was calculated
    :type freq: float
    :return: Full path to written polyco file
    :rtype: str
    """
    with open("craftfrb.polyco", "w") as polycoout:
        polycoout.write(
            f"fake+fake DD-MMM-YY %02d%02d%05.2f %.15f %.4f 0.0 0.0\n"
            % (hh, mm, ss, polycorefmjd, dm)
        )
        polycoout.write(
            f"0.0 {1.0/float(fakepulsarperiod):.3f} 0 100 3 {freq:.3f}\n"
        )
        polycoout.write(
            "0.00000000000000000E-99 "
            "0.00000000000000000E-99 "
            "0.00000000000000000E-99\n"
        )
        polycoout.close()

    return f"{os.curdir}/craftfrb.polyco"


def write_binconfig(
    fname: str,
    polycopath: str,
    binedges: "list[float]",
    weights: "list[float]",
    scrunch: bool = True,
) -> None:
    """Write a binconfig file.

    :param fname: File name to save binconfig in
    :type fname: str
    :param polycopath: Path to the polyco file associated with the bin
        config
    :param binedges: Phase ends of each bin (i.e. bin edges in units
        of pulse phase)
    :type phaseends: list[float]
    :param scrunch: Set the SCRUNCH OUTPUT parameter of the bin config
        to True or False [Default = True]
    :type scrunch: bool, optional
    :param weights: Bin weights as floats between 0 and 1.
    :type weights: list[float]
    """
    assert len(binedges) == len(
        weights
    ), f"{fname}: Must provide same number of bin edges as weights"

    nbins = len(binedges)

    with open(fname, "w") as binconfout:
        binconfout.write(f"NUM POLYCO FILES:   1\n")
        binconfout.write(f"POLYCO FILE 0:      {polycopath}\n")
        binconfout.write(f"NUM PULSAR BINS:    {nbins}\n")

        if scrunch:
            binconfout.write(f"SCRUNCH OUTPUT:     TRUE\n")
        else:
            binconfout.write(f"SCRUNCH OUTPUT:     FALSE\n")

        for i in range(nbins):
            binconfout.write(f"BIN PHASE END {i}:    {binedges[i]:.9f}\n")
            binconfout.write(f"BIN WEIGHT {i}:       {weights[i]}\n")


def calc_gate_bins(
    cand: "list[str]",
    timediff: float,
    polycorefmjd: float,
) -> "tuple[list[float], list[float]]":
    """Determine the bins for the gate binconfig

    The gate mode has two bins:
        > On-pulse (from 50 milliseconds before to 150 milliseconds
            after the burst)
        > Off-pulse (everything else) - weighted 0

    :param cand: Fields of the snoopy candidate
    :type cand: list[str]
    :param timediff: The time difference between the VCRAFT and snoopy
        log arrival times for the pulse, including geometric delay, in
        s
    :type timediff: float
    :param polycorefmjd: Polyco reference time in MJD
    :type polycorefmjd: float
    :return: List of bin edges and list of weights
    :rtype: tuple[list[float], list[float]]
    """
    pulsewidthms = float(cand[3]) * 1.7
    mjd = float(cand[7])

    gatestartmjd = mjd - (pulsewidthms + 100) / (
        2 * 86400000.0
    )  # pulse width is in ms at this point
    gateendmjd = (
        gatestartmjd + (pulsewidthms + 200) / 86400000.0
    )  # pulse width is in ms at this point
    gatestartphase = (
        86400.0 * (gatestartmjd - polycorefmjd) + timediff
    ) / fakepulsarperiod
    gateendphase = (
        86400.0 * (gateendmjd - polycorefmjd) + timediff
    ) / fakepulsarperiod

    binedges = [gatestartphase, gateendphase]
    weights = [0.0, 1.0]

    return binedges, weights


def calc_rfi_bins(
    gatebinedges: "list[float]",
) -> "tuple[list[float], list[float]]":
    """Determine the bins for the RFI binconfig

    The RFI mode creates a bin before and after the expected position of
    the FRB based on the bin edges from the gate mode. This is to give
    data that should contain only RFI that can then be subtracted from
    the on-signal data to give relatively RFI-free data.

    :param gatebinedges: Bin edges of the gate binconfig
    :type gatebinedges: list[float]
    :return: List of bin edges and list of weights
    :rtype: tuple[list[float], list[float]]
    """
    gatestartphase = gatebinedges[0]
    gateendphase = gatebinedges[1]

    rfistartphase1 = (
        gatestartphase - 0.02 / fakepulsarperiod
    )  # RFI gate (early side) starts 20ms before the start of the pulse
    rfiendphase1 = (
        gatestartphase - 0.004 / fakepulsarperiod
    )  # RFI gate (early side) ends 4ms before the start of the pulse
    rfistartphase2 = (
        gateendphase + 0.004 / fakepulsarperiod
    )  # RFI gate (late side) starts 4ms after the end of the pulse
    rfiendphase2 = (
        gateendphase + 0.02 / fakepulsarperiod
    )  # RFI gate (early side) ends 20ms after the end of the pulse

    gateedges = [rfistartphase1, rfiendphase1, rfistartphase2, rfiendphase2]
    binedges = [0.0, 1.0, 0.0, 1.0]

    return gateedges, binedges


def calc_htr_bins(
    cand: "list[str]",
    gatestartphase: float,
) -> "tuple[list[float], list[float]]":
    """Determine the bins for the high time resolution binconfig

    The high time resolution mode creates bins 216 us wide in the range
    2 ms before and after the expected FRB position (as determined by
    the gate bin edges).

    :param cand: Fields of the snoopy candidate
    :type cand: list[str]
    :param gatestartphase: Start phase of the gate bin
    :type gatestartphase: float
    :return: List of bin edges and list of weights
    :rtype: tuple[list[float], list[float]]
    """
    pulsewidthms = float(cand[3]) * 1.7
    mjd = float(cand[7])

    gatestartmjd = mjd - (pulsewidthms + 1000) / (
        2 * 86400000.0
    )  # pulse width is in ms at this point

    binstartmjd = mjd - fakepulsarperiod * pulsewidthms / (
        2 * 86400000.0
    )  # pulse width is in ms at this point
    gateendmjd = (
        gatestartmjd + pulsewidthms / 86400000.0
    )  # pulse width is in ms at this point
    binmicrosec = 216
    extrawidth = 2  # ms on either side of the snoopy detected pulse
    binstartphase = gatestartphase - float(extrawidth) / (
        1000.0 * fakepulsarperiod
    )
    bindeltaphase = binmicrosec / (fakepulsarperiod * 1e6)
    numbins = int((pulsewidthms + 2 * extrawidth) / (binmicrosec / 1000.0))

    binedges = [binstartphase + i * bindeltaphase for i in range(numbins + 1)]
    binweights = [1 for i in range(numbins + 1)]

    return binedges, binweights, numbins


def calc_finder_bins(
    rfibinedges: "list[float]", numfinderbins: int = 20
) -> "tuple[list[float], list[float]]":
    """Determine the bins for the finder binconfig

    The finder mode creates (by default) 20 bins over the same range as
    the single on-pulse bin in the RFI mode. These will normally be
    around 10 milliseconds wide each.

    :param rfibinedges: Bin edges of the rfi binconfig
    :type rfibinedges: list[float]
    :param numfinderbins: Number of finder bins [Default = 20]
    :type numfinderbins: int, optional
    :return: List of bin edges and list of weights
    :rtype: tuple[list[float], list[float]]
    """
    rfistartphase1 = rfibinedges[0]
    rfiendphase1 = rfibinedges[1]
    rfistartphase2 = rfibinedges[2]
    rfiendphase2 = rfibinedges[3]

    bindeltaphase = (rfistartphase2 - rfiendphase1) / numfinderbins
    binstartphase = rfiendphase1

    binedges = [
        binstartphase + i * bindeltaphase for i in range(numfinderbins + 1)
    ]
    binweights = [1 for i in range(numfinderbins + 1)]

    return binedges, binweights, numfinderbins


def write_subtractions_script(
    gatebinedges: "list[float]",
    rfibinedges: "list[float]",
    finderbinedges: "list[float]",
    htrbinedges: "list[float]",
    numbins: int,
    numfinderbins: int,
) -> None:
    """Write the bash script that performs the RFI subtractions.

    TODO: understand this function - what are the scales?

    :param gatebinedges: Bin edges for the gate mode
    :type gatebinedges: list[float]
    :param rfibinedges: Bin edges for the RFI mode
    :type rfibinedges: list[float]
    :param finderbinedges: Bin edges for the finder mode
    :type finderbinedges: list[float]
    :param htrbinedges: Bin edges for the high time resolution mode
    :type htrbinedges: list[float]
    :param numbins: Number of high time resolution bins
    :type numbins: int
    :param numfinderbins: Number of finder bins
    :type numfinderbins: int
    """
    binscale = (htrbinedges[1] - htrbinedges[0]) / (
        rfibinedges[3] + rfibinedges[1] - rfibinedges[2] - rfibinedges[0]
    )

    finderbinscale = (finderbinedges[1] - finderbinedges[0]) / (
        rfibinedges[3] + rfibinedges[1] - rfibinedges[2] - rfibinedges[0]
    )

    gatescale = (gatebinedges[1] - gatebinedges[0]) / (
        rfibinedges[3] + rfibinedges[1] - rfibinedges[2] - rfibinedges[0]
    )
    with open("dosubtractions.sh", "w") as subout:
        subout.write(
            f"uvsubScaled.py *_gate.fits *_rfi.fits {gatescale:.9f} gate_norfi.fits\n"
        )
        for i in range(numbins):
            subout.write(
                f"uvsubScaled.py *_bin{i:02d}.fits *_rfi.fits {binscale:.9f} bin{i:02d}_norfi.fits\n"
            )
        for i in range(numfinderbins):
            subout.write(
                f"uvsubScaled.py finderbin{i:02d}.fits *_rfi.fits {finderbinscale:.9f} finderbin{i:02d}_norfi.fits\n"
            )
        subout.close()


if __name__ == "__main__":
    _main()
