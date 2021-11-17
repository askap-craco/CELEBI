#!/usr/bin/env ParselTongue
# Imports ########################################################
import argparse
import datetime
import glob
import os
import socket
import sys

import astropy.units as un
import vlbatasks
from AIPS import AIPS
from AIPSData import AIPSImage, AIPSUVData
from AIPSTask import AIPSTask
from astropy.time import Time

# Global variables
try:
    AIPSVER = os.environ["PSRVLBAIPSVER"]
except KeyError:
    AIPSVER = "31DEC18"

AIPSDISK = 1
DELAYWINDOW = 0  # Search everything
RATEWINDOW = -1  # Don't search rates
PLOTSPERPAGE = 4
OUTKLASS = "SPLIT"
SOLINTMINS = 1  # Long enough that we just get one solutions
SUMIFS = False
SEQNO = 1


def _main():
    args = get_args()

    AIPS.userno = args.userno
    xpolmodelfile = args.xpoldelaymodelfile

    # TODO: reformat with these as global variables
    snversion = 1
    clversion = 1
    bpversion = 1

    if xpolmodelfile != "":
        xpol_prefix = "_xpol"
    else:
        xpol_prefix = "_noxpol"
    if args.src != "":
        src = "_" + args.src
    else:
        src = args.src

    # Make path names absolute if needed
    targetpath = os.path.abspath(args.target)
    calpath = os.path.abspath(args.calibrator)
    if xpolmodelfile != "":
        xpolmodelfile = os.path.abspath(xpolmodelfile)

    # Calibrated output filenames
    targetoutfname, targetmsfname = out_fnames(targetpath)
    caloutfname, calmsfname = out_fnames(calpath)

    # Solution output filenames
    calsolnfnames = soln_fnames(xpolmodelfile, args.src)

    # If we are running targetonly, check that all the calibration files exist
    if args.targetonly:
        validate_soln_files(calsolnfnames)

    bpfname = calsolnfnames[0]
    fringsnfname = calsolnfnames[1]
    selfcalsnfname = calsolnfnames[2]
    xpolsnfname = calsolnfnames[3]
    bptableplotfname = calsolnfnames[4]
    uncalxcorplotfname = calsolnfnames[5]
    allcalxcorplotfname = calsolnfnames[6]

    # Load and flag the target data if needed
    if not args.calibrateonly:
        targetdata = load_data(args.target, args.uvsrt)
        if args.tarflagfile != "":
            flag_data(targetdata, args.tarflagfile, args.shadow)

    # Load and flag the calibrator data if needed
    if not args.targetonly:
        caldata = load_data(args.calibrator, args.uvsrt)
        if args.flagfile != "":
            flag_data(caldata, args.flagfile, args.shadow)

        # Get the reference frequency of the dataset
        reffreqs = get_ref_freqs(caldata)

    # Run CLCOR to correct PANG if needed
    if xpolmodelfile != "":
        if not args.targetonly:
            vlbatasks.clcor_pang(caldata, clversion)
        if not args.calibrateonly:
            vlbatasks.clcor_pang(targetdata, clversion)
        clversion = clversion + 1

    # Run FRING
    if not args.targetonly:
        run_FRING(
            caldata,
            snversion,
            clversion,
            args.sourcename,
            args.refant,
            fringsnfname,
        )

    # Load FRING SN table into the target
    if not args.calibrateonly:
        vlbatasks.loadtable(targetdata, fringsnfname, snversion)

    # Calibrate
    if not args.targetonly:
        vlbatasks.applysntable(
            caldata, snversion, "SELN", clversion, args.refant
        )
    if not args.calibrateonly:
        vlbatasks.applysntable(
            targetdata, snversion, "SELN", clversion, args.refant
        )

    snversion += 1
    clversion += 1

    # Correct for leakage if needed
    if xpolmodelfile != "":
        # First the xpoldelays
        if not args.targetonly:
            correct_leakage(
                caldata,
                snversion,
                clversion,
                args.sourcename,
                args.refant,
                xpolmodelfile,
                xpolsnfname,
            )
        if not args.calibrateonly:
            vlbatasks.loadtable(targetdata, xpolsnfname, snversion)
            vlbatasks.applysntable(
                targetdata, snversion, "2PT", clversion, args.refant
            )
        snversion += 1
        clversion += 1

    # Run bandpass correction
    if not args.targetonly:
        run_bandpass(
            caldata,
            clversion,
            bpversion,
            args.sourcename,
            bpfname,
            args.cpasspoly,
            args.bpass,
        )
        # Plot the bandpass table
        if not args.skipplot:
            plot_bandpass(caldata, bpversion, bptableplotfname)

    # Load up the bandpass to the target
    if not args.calibrateonly:
        vlbatasks.loadtable(targetdata, bpfname, bpversion)

    # Run selfcal
    if not args.targetonly:
        run_selfcal(
            caldata,
            clversion,
            args.sourcename,
            args.refant,
            args.flux,
            selfcalsnfname,
        )

    # Load up the selfcal SN table
    if not args.calibrateonly:
        vlbatasks.loadtable(targetdata, selfcalsnfname, snversion)
    if not args.targetonly:
        vlbatasks.loadtable(caldata, selfcalsnfname, snversion)

    # Calibrate
    if not args.targetonly:
        vlbatasks.applysntable(
            caldata, snversion, "SELN", clversion, args.refant
        )
    if not args.calibrateonly:
        vlbatasks.applysntable(
            targetdata, snversion, "SELN", clversion, args.refant
        )
    snversion += 1
    clversion += 1

    # Plot the uncalibrated and calibrated cross-correlation results if desired
    if not args.targetonly:
        if not args.skipplot:
            plot_xcor(
                caldata,
                clversion,
                bpversion,
                args.xcorplotsmooth,
                uncalxcorplotfname,
                allcalxcorplotfname,
            )

    # Run SPLIT and write output data for calibrator
    if not args.targetonly:
        run_split(caldata, caloutfname, clversion, args.sourcename)

    # Run SPLIT and write output data for target
    if not args.calibrateonly:
        run_split(targetdata, targetoutfname, clversion, args.sourcename)

    # Create a README file for the calibration and a tarball with it plus all the calibration
    if not args.targetonly:
        write_readme(calsolnfnames, args.calibrator, xpolmodelfile, reffreqs)

    # Convert to a measurement set
    if not args.targetonly:
        casaout = open("loadtarget.py", "w")
        casaout.write(
            f"importuvfits(fitsfile='{caloutfname}',vis='{calmsfname}',antnamescheme='old')\n"
        )
        casaout.close()
        os.system("runloadtarget.sh")

    if not args.calibrateonly:
        casaout = open("loadtarget.py", "w")
        casaout.write(
            f"importuvfits(fitsfile='{targetoutfname}',vis='{targetmsfname}',antnamescheme='old')\n"
        )
        casaout.close()
        os.system("runloadtarget.sh")

    def write_timestep_loop(write_file, fmt_str, vals):
        # iterate over the duration given in the options, writing the fmt_str to
        # write_file (formatted with vals), with the addition of selecting the data
        # for the step and appending an index to the imagename
        init_imagename = vals[1]

        i = 0
        start = Time(args.start, format="mjd")
        while (i + 1) * args.inttime < args.duration:
            # get end time by adding inttime to start time
            end = start + args.inttime * un.s

            # format times as YYYY/MM/DD/hh:mm:ss strings
            start_str = start.iso.replace(" ", "/").replace("-", "/")
            end_str = end.iso.replace(" ", "/").replace("-", "/")

            # add data selection to fmt_str
            new_fmt_str = fmt_str.replace(
                ")", f", selectdata=True, timerange='{start_str}~{end_str}')\n"
            )

            # specify image name for this timestep
            vals[1] = init_imagename + "_" + str(i)

            write_file.write(new_fmt_str.format(*vals))

            i += 1
            start = end

        print("Wrote " + str(i) + " time steps")

    def write_single_timestep(write_file, fmt_str, vals, i):
        start = Time(args.start, format="mjd") + args.inttime * un.s * i
        end = start + args.inttime * un.s
        start_str = start.iso.replace(" ", "/").replace("-", "/")
        end_str = end.iso.replace(" ", "/").replace("-", "/")

        # add data selection to fmt_str
        new_fmt_str = fmt_str.replace(
            ")", f", selectdata=True, timerange='{start_str}~{end_str}')\n"
        )

        # specify image name for this timestep
        vals[1] = f"{i:03d}_{vals[1]}"

        write_file.write(new_fmt_str.format(*vals))

    # Run the imaging via CASA if desired
    if not args.calibrateonly:
        do_imaging(args)

    imagesize = args.imagesize
    polarisations = args.pols.split(",")

    if args.dirtyonly:
        polarisations = ["I"]

    if args.dirtymfs:
        polarisations = ["I"]

    if args.imagecube:
        # Do the cube
        for pol in polarisations:
            casaout = open("imagescript.py", "w")
            imagename = f"TARGET.cube.{pol}"
            offsourcename = f"OFFSOURCE.cube.{pol}"
            maskstr = "'circle[[{0}pix,{0}pix] ,5pix ]'".format(imagesize / 2)
            phasecenter = f"'{args.phasecenter}'".encode()
            imsize = "[{0},{0}]".format(imagesize)

            # If desired, produce the noise image
            if args.noisecentre:
                rmscenter = f"{args.noisecentre}"
                rmsimsize = "[{0},{0}]".format(imagesize * 4)
                os.system(f"rm -rf {offsourcename}*")
                os.system(f"rm -rf {imagename}*")
                outlierfile = open(f"outlierfield_Stokes{pol}.txt", "w")
                outlierfile.write(
                    "imagename={0}\nimsize={1}\nphasecenter={2}\nmask="
                    "circle[[{3}pix,{3}pix] ,{3}pix ]"
                    "\n".format(
                        offsourcename, rmsimsize, rmscenter, imagesize * 2
                    )
                )
                outlierfile.close()
                outlierfields = f"'outlierfield_Stokes{pol}.txt'"
            else:
                outlierfields = "[]"
                os.system(f"rm -rf {imagename}.*")

            # If desired, produce the only the dirty image
            if args.dirtyonly:
                imagename = f"TARGET.cube.dirim.{pol}"
                os.system(f"rm -rf {imagename}.*")
                fmt_str = "tclean(vis='{0}', imagename='{1}', imsize={2}, cell=['{8}arcsec', '{8}arcsec'], stokes='{3}', specmode='cube', width={4}, phasecenter={5}, gridder='widefield', wprojplanes=-1, pblimit=-1, deconvolver='multiscale', weighting='natural', niter=0, mask={6}, outlierfile={7})"
                vals = [
                    targetmsfilename,
                    imagename,
                    imsize,
                    pol,
                    args.averagechannels,
                    phasecenter,
                    maskstr,
                    outlierfields,
                    pixelsize,
                ]
                if args.doalltimesteps:
                    write_timestep_loop(casaout, fmt_str, vals)
                elif args.dosingletimestep:
                    write_single_timestep(
                        casaout, fmt_str, vals, args.timestep
                    )
                else:
                    casaout.write(
                        "tclean(vis='{0}', imagename='{1}', imsize={2}, cell=['{8}arcsec', '{8}arcsec'], stokes='{3}', specmode='cube', width={4}, phasecenter={5}, gridder='widefield', wprojplanes=-1, pblimit=-1, deconvolver='multiscale', weighting='natural', niter=0, mask={6}, outlierfile={7})".format(
                            targetmsfilename,
                            imagename,
                            imsize,
                            pol,
                            args.averagechannels,
                            phasecenter,
                            maskstr,
                            outlierfields,
                            pixelsize,
                        )
                    )

            elif args.dirtymfs:
                imagename = f"TARGET.mfs.dirim.{pol}"
                os.system(f"rm -rf {imagename}.*")
                fmt_str = "tclean(vis='{0}', imagename='{1}', imsize={2}, cell=['{8}arcsec', '{8}arcsec'], stokes='{3}', specmode='mfs', width={4}, phasecenter={5}, gridder='widefield', wprojplanes=-1, pblimit=-1, deconvolver='multiscale', weighting='natural', niter=0, mask={6}, outlierfile={7})"
                vals = [
                    targetmsfilename,
                    imagename,
                    imsize,
                    pol,
                    args.averagechannels,
                    phasecenter,
                    maskstr,
                    outlierfields,
                    pixelsize,
                ]
                if args.doalltimesteps:
                    write_timestep_loop(casaout, fmt_str, vals)
                elif args.dosingletimestep:
                    write_single_timestep(
                        casaout, fmt_str, vals, args.timestep
                    )

                casaout.write(
                    "tclean(vis='{0}', imagename='{1}', imsize={2}, cell=['{8}arcsec', '{8}arcsec'], stokes='{3}', specmode='mfs', width={4}, phasecenter={5}, gridder='widefield', wprojplanes=-1, pblimit=-1, deconvolver='multiscale', weighting='natural', niter=0, mask={6}, outlierfile={7})".format(
                        targetmsfilename,
                        imagename,
                        imsize,
                        pol,
                        args.averagechannels,
                        phasecenter,
                        maskstr,
                        outlierfields,
                        pixelsize,
                    )
                )

            elif args.cleanmfs:
                imagename = args.imagename
                os.system(f"rm -rf {imagename}.*")
                fmt_str = "tclean(vis='{0}', imagename='{1}', imsize={2}, cell=['{7}arcsec', '{7}arcsec'], stokes='{3}', specmode='mfs', phasecenter={4}, gridder='widefield', wprojplanes=-1, pblimit=-1, deconvolver='multiscale', weighting='natural', niter=100, mask={5}, outlierfile={6}, savemodel='modelcolumn')"
                vals = [
                    targetmsfilename,
                    imagename,
                    imsize,
                    pol,
                    phasecenter,
                    maskstr,
                    outlierfields,
                    pixelsize,
                ]
                if args.doalltimesteps:
                    write_timestep_loop(casaout, fmt_str, vals)
                elif args.dosingletimestep:
                    write_single_timestep(
                        casaout, fmt_str, vals, args.timestep
                    )
                else:
                    casaout.write(
                        "tclean(vis='{0}', imagename='{1}', imsize={2}, cell=['{7}arcsec', '{7}arcsec'], stokes='{3}', specmode='mfs', phasecenter={4}, gridder='widefield', wprojplanes=-1, pblimit=-1, deconvolver='multiscale', weighting='natural', niter=100, mask={5}, outlierfile={6}, savemodel='modelcolumn')".format(
                            targetmsfilename,
                            imagename,
                            imsize,
                            pol,
                            phasecenter,
                            maskstr,
                            outlierfields,
                            pixelsize,
                        )
                    )

            # Default: produce a cleaned cube image
            else:
                fmt_str = "tclean(vis='{0}', imagename='{1}', imsize={2}, cell=['{8}arcsec', '{8}arcsec'], stokes='{3}', specmode='cube', width={4}, gridder='widefield', wprojplanes=-1, pblimit=-1, deconvolver='multiscale', weighting='natural', niter=5000, cycleniter=100, mask={5}, savemodel='modelcolumn', phasecenter={6}, outlierfile={7}, spw={9})"
                vals = [
                    targetmsfilename,
                    imagename,
                    imsize,
                    pol,
                    args.averagechannels,
                    maskstr,
                    phasecenter,
                    outlierfields,
                    pixelsize,
                    args.spwrange,
                ]
                if args.doalltimesteps:
                    write_timestep_loop(casaout, fmt_str, vals)
                elif args.dosingletimestep:
                    write_single_timestep(
                        casaout, fmt_str, vals, args.timestep
                    )
                else:
                    casaout.write(
                        "tclean(vis='{0}', imagename='{1}', imsize={2}, cell=['{8}arcsec', '{8}arcsec'], stokes='{3}', specmode='cube', width={4}, gridder='widefield', wprojplanes=-1, pblimit=-1, deconvolver='multiscale', weighting='natural', niter=5000, cycleniter=100, mask={5}, savemodel='modelcolumn', phasecenter={6}, outlierfile={7}, spw={9})".format(
                            targetmsfilename,
                            imagename,
                            imsize,
                            pol,
                            args.averagechannels,
                            maskstr,
                            phasecenter,
                            outlierfields,
                            pixelsize,
                            args.spwrange,
                        )
                    )

            casaout.close()
            os.system("chmod 775 imagescript.py")
            os.system("runimagescript.sh")

            # If desired, also export the image as a FITS file
            if args.exportfits:
                casaout = open("exportfits.py", "w")
                casaout.write(
                    'exportfits(imagename="{0}.image",fitsimage="{0}.FITS")\n'.format(
                        imagename
                    )
                )
                casaout.close()
                os.system("chmod 775 exportfits.py")
                os.system("runexportfits.sh")

            # If desired, also make the JMFIT output
            if args.imagejmfit:
                casaout = open("imagescript.py", "w")
                casaout.write(
                    f'exportfits(imagename="{imagebase}.image",fitsimage="{imagebase}.fits")\n'
                )
                casaout.close()
                os.system(
                    "casa --nologger -c imagescript.py"
                )  # Get the number of channels in the dataset
                numchannels = vlbatasks.getNumChannels(targetdata)
                for i in range(numchannels / args.averagechannels):
                    locstring = "%d,%d,%d,%d,%d,%d" % (
                        imagesize / 2 - 12,
                        imagesize / 2 - 12,
                        i,
                        imagesize / 2 + 12,
                        imagesize / 2 + 12,
                        i,
                    )
                    os.system(
                        "jmfitfromfile.py %s.fits %s.slice%03d.jmfit.stats %s"
                        % (imagebase, imagebase, i, locstring)
                    )


def get_args() -> argparse.Namespace:
    """Parse command line arguments

    :return: Command line argument paramters
    :rtype: :class:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_option(
        "-t",
        "--target",
        default="",
        help="The target FITS file, to be calibrated and optionally imaged",
    )
    parser.add_option(
        "-c",
        "--calibrator",
        default="",
        help="The calibrator FITS file (usually on 0407)",
    )
    parser.add_option(
        "-r",
        "--refant",
        type=int,
        default=3,
        help="The reference antenna to use.",
    )
    parser.add_option(
        "-u",
        "--userno",
        type=int,
        default=2,
        help="The AIPS user number",
    )
    parser.add_option(
        "-s",
        "--sourcename",
        default="CRAFTSRC",
        help="The name of the source in the FITS files",
    )
    parser.add_option(
        "-f",
        "--flux",
        type=float,
        default=9.5,  # 0407 flux
        help="Calibrator flux in Jy,  Defaulted to correct value for 0407",
    )
    parser.add_option(
        "-i",
        "--imagecube",
        default=False,
        action="store_true",
        help="Image the FRB in chunks (i.e., make a cube)",
    )
    parser.add_option(
        "--dirtyonly",
        default=False,
        action="store_true",
        help="Just make a dirty image, cube",
    )
    parser.add_option(
        "--dirtymfs",
        default=False,
        action="store_true",
        help="Just make a dirty image, mfs",
    )
    parser.add_option(
        "--cleanmfs",
        default=False,
        action="store_true",
        help="Make a clean image, mfs",
    )
    parser.add_option(
        "--imagename",
        default="TARGET",
        help="The name of the image files for the cleaned MFS images. Default:"
        " TARGET",
    )
    parser.add_option(
        "--exportfits",
        default=False,
        action="store_true",
        help="Export image as FITS file",
    )
    parser.add_option(
        "--spwrange",
        type=str,
        default="''",
        help="The spw range used for imaging if selected. CASA format: '0:0~42' for spectral window 0, channels 0 to 42. Default is all channels.",
    )
    parser.add_option(
        "--calibrateonly",
        default=False,
        action="store_true",
        help="Only generate the calibration files, don't do anything with target",
    )
    parser.add_option(
        "--targetonly",
        default=False,
        action="store_true",
        help="Use saved calibration files",
    )
    parser.add_option(
        "--bpass",
        default=False,
        action="store_true",
        help="Use BPASS rather than CPASS to do the bandpass correction",
    )
    parser.add_option(
        "-j",
        "--imagejmfit",
        default=False,
        action="store_true",
        help="Jmfit the individual slices of the cube",
    )
    parser.add_option(
        "--cpasspoly",
        default=10,
        type=int,
        help="Number of polynomial terms in CPASS",
    )
    parser.add_option(
        "-a",
        "--averagechannels",
        type=int,
        default=24,
        help="Number of channels to average together per cube slice",
    )
    parser.add_option(
        "-F",
        "--flagfile",
        default="",
        help="Flag file to apply to calibrator data only, if desired. Used to ensure RFI doesn't corrupt FRING or BPASS.",
    )
    parser.add_option(
        "-g",
        "--tarflagfile",
        default="",
        help="Flag file to apply to target data only, if desired. Used to flag any necessary channels for, e.g., RFI or missing data",
    )
    parser.add_option(
        "--shadow",
        nargs=2,
        default=None,
        help="Set if flagging due to shadowing is desired. Takes two arguments: arg1 > 0 flag for shadowing; shadow diameter in m. arg2 flag for cross-talk; baseline (BL) in m",
    )
    parser.add_option(
        "-p",
        "--phasecenter",
        default="",
        help="phase center for the target field (blank will leave it at correlation centre)",
    )
    # parser.add_option("-l", "--leakagecorrect", default=False, action="store_true",
    #                  help="Run lpcal to try and correct any leakage present")
    parser.add_option(
        "-x",
        "--xpoldelaymodelfile",
        default="",
        help="Model to use for xpol delay correction (blank = no correction)",
    )
    parser.add_option(
        "--imagesize", type=int, default=128, help="Size of the image to make"
    )
    parser.add_option(
        "--xcorplotsmooth",
        type=int,
        default=32,
        help="Length of the smoothing kernel in channels for xcor plotting",
    )
    parser.add_option(
        "--skipplot",
        default=False,
        action="store_true",
        help="Skip the plotting to save time",
    )
    parser.add_option(
        "--pixelsize", type=float, default=1, help="Pixel size in arcseconds"
    )
    parser.add_option(
        "--uvsrt",
        default=False,
        action="store_true",
        help="Run UVSRT on the data after loading",
    )
    parser.add_option(
        "--noisecentre",
        default="",
        help="CASA format position at which noise should be estimated, blank=don't make an off-source image",
    )
    parser.add_option(
        "--src", default="", help="Name of the target (e.g., FRB or Vela)"
    )
    parser.add_option(
        "--pols",
        type=str,
        default="XX,YY,I,Q,U,V",
        help='The polarisations to be imaged if --imagecube is set. Defaulted to all. Input as a list of strings: e.g., "XX,YY"',
    )
    parser.add_option(
        "--doalltimesteps",
        default=False,
        action="store_true",
        help="Image all time steps separately",
    )
    parser.add_option(
        "--dosingletimestep",
        default=False,
        action="store_true",
        help="Image a single time step",
    )
    parser.add_option(
        "--timestep",
        type=int,
        default=None,
        help="Index of time step to image",
    )
    parser.add_option(
        "--start", type=float, default=None, help="Start of time range as MJD"
    )
    parser.add_option(
        "--inttime", type=float, default=None, help="Time step duration in s"
    )
    parser.add_option(
        "--duration",
        type=float,
        default=None,
        help="Duration over which to iterate time steps",
    )

    return parser.parse_args()


def out_fnames(fitspath: str) -> "tuple[str, str]":
    """Determine output filename based on the provided path to a target
    or calibrator fits file.

    Replaces the suffix (either `.fits` or `.uvfits`) with
    `_calibrated_uv.fits`, and the path to the provided file with the
    current working directory's path. Also returns a path with the same
    structure, but with `.ms` as the extension instead of `.fits.`

    Additionally checks if the `.ms` file already exists. If it does,
    abort the program.

    :param fitspath: Path to the input fits file. Must have either
        `.fits` or `.uvfits` as its suffix
    :type fitspath: str
    :return: Path for two output files: a `.fits` file and a `.ms`
        measurement set.
    :rtype: str
    """
    assert (
        fitspath[:-5] == ".fits" or fitspath[:-7] == ".uvfits"
    ), "fits filename must end in either .fits or .uvfits!"

    fitsbase = fitspath.split("/")[-1]
    if ".uvfits" in fitspath:
        fitsfname = f"{os.getcwd}/{fitsbase[:-7]}_calibrated_uv.fits"
    else:
        fitsfname = f"{os.getcwd}/{fitsbase[:-5]}_calibrated_uv.fits"

    msfname = fitsfname[:-4] + "ms"

    if os.path.exists(msfname):
        print(f"{msfname} already exists - aborting!!!")
        sys.exit()

    return fitsfname, msfname


def soln_fnames(xpol_prefix: str, src: str) -> "tuple[str]":
    """Determine filenames for the calibration solutions.

    :param xpol_prefix: `_xpol` if we have an xpolmodelfile, `_noxpol`
        otherwise
    :type xpol_prefix: str
    :param src: Source name (i.e. FRB/Vela/etc.) with a "_" pre-appended
    :type src: str
    :return: A prefix specifilenames for:
            - Bandpasses
            - Fring delays
            - Selcal solutions
            - X polarisation fring delays
            - Bandpass table plots
            - Uncalibrated cross-correlation plots
            - Calibrated cross-correlation plots
            - README
            - Tarball for calibration solutions
    :rtype: tuple[str]
    """
    bpfname = os.path.abspath(f"bandpasses{xpol_prefix}{src}.bp.txt")
    fringsnfname = os.path.abspath(f"delays{xpol_prefix}{src}.sn.txt")
    selfcalsnfname = os.path.abspath(f"selfcal{xpol_prefix}{src}.sn.txt")
    xpolsnfname = os.path.abspath(f"xpolfring{xpol_prefix}{src}.sn.txt")
    bptableplotfname = os.path.abspath(f"bptable{xpol_prefix}{src}.ps")
    uncalxcorplotfname = os.path.abspath(f"uncalxcor{xpol_prefix}{src}.ps")
    allcalxcorplotfname = os.path.abspath(f"allcalxcor{xpol_prefix}{src}.ps")
    readmefname = os.path.abspath(f"README{xpol_prefix}{src}.calibration")
    calibtarballfile = f"calibration{xpol_prefix}{src}.tar.gz"

    return (
        xpol_prefix,
        src,
        bpfname,
        fringsnfname,
        selfcalsnfname,
        xpolsnfname,
        bptableplotfname,
        uncalxcorplotfname,
        allcalxcorplotfname,
        readmefname,
        calibtarballfile,
    )


def validate_soln_files(
    calsolnfnames: "tuple[str, str, str, str, str, str, str]",
) -> None:
    """Check that the calibration solution files already exist. If
    any mission-critical files are missing (i.e. not including plot
    files), abort.

    This should only be done when running the `targetonly` mode.

    :param calsolnfnames: Filenames for:
        - Bandpasses
        - Fring delays
        - Selcal solutions
        - X polarisation fring delays
        - Bandpass table plots
        - Uncalibrated cross-correlation plots
        - Calibrated cross-correlation plots
    :type cal_soln_fnames: tuple[str, str, str, str, str, str, str]
    """
    bpfname = calsolnfnames[0]
    fringsnfname = calsolnfnames[1]
    selfcalsnfname = calsolnfnames[2]
    xpolsnfname = calsolnfnames[3]
    missingfiles = []
    if not os.path.exists(bpfname):
        missingfiles.append(bpfname)
    if not os.path.exists(fringsnfname):
        missingfiles.append(fringsnfname)
    if not os.path.exists(selfcalsnfname):
        missingfiles.append(selfcalsnfname)
    if not os.path.exists(xpolsnfname):
        missingfiles.append(xpolsnfname)
    if len(missingfiles) > 0:
        print(
            "Running targetonly but the following files are missing:",
            missingfiles,
        )
        sys.exit()


def load_data(fits: str, uvsrt: bool):
    """Load data from a fits file using AIPS.

    TODO: Figure out what type data is returned as!

    :param fits: Filename of the visibility fits file
    :type fits: str
    :param uvsrt: If True, runs UVSRT on the data after loading
    :type uvsrt: bool
    """
    data = vlbatasks.zapAndCreateUVData("CRAFTTARG", "UVDATA", AIPSDISK, 1)
    vlbatasks.fitld_corr(fits, data, [], "", 0.0001)

    if uvsrt:
        sorteddata = vlbatasks.zapAndCreateUVData(
            "CRAFTTARG", "UVSRT", AIPSDISK, 1
        )
        vlbatasks.uvsrt(data, sorteddata)
        data.zap()
        data = sorteddata

    return data


def flag_data(data, flagfile: str, shadow: "list[str]") -> None:
    """Flag data according to provided flagfile

    :param data: Data to be flagged
    :type data: [type]
    :param flagfile: Path to file containing flags to be applied
    :type flagfile: str
    :param shadow: Shadow command line argument. If not None, then apply
        shadowing flags
    :type shadow: list[str]
    :return: Flagged data
    :rtype: [type]
    """
    vlbatasks.userflag(data, 1, flagfile)

    if shadow:
        shadowdiameter = float(shadow[0])
        xtalkbl = float(shadow[1])
        vlbatasks.shadowflag(data, 1, shadowdiameter, xtalkbl)
        print("Shadowing diameter: " + str(shadowdiameter))
        print("Cross-talk baseline: " + str(xtalkbl))


def get_ref_freqs(caldata) -> "list[float]":
    """Get the reference frequencies from the calibrator data

    :param caldata: Calibrator data
    :type caldata: [type]
    :return: List of reference frequencies
    :rtype: list[float]
    """
    reffreqs = []
    fqtable = caldata.table("FQ", 1)
    for row in fqtable:
        try:
            for iffreq in row.if_freq:
                freqentry = float(iffreq) + float(caldata.header.crval[2])
                reffreqs.append(float(iffreq) + float(caldata.header.crval[2]))
        except (AttributeError, TypeError):
            freqentry = float(row.if_freq) + float(caldata.header.crval[2])
            reffreqs.append(
                float(row.if_freq) + float(caldata.header.crval[2])
            )
    return reffreqs


def run_FRING(
    caldata,
    snversion: int,
    clversion: int,
    sourcename: str,
    refant: int,
    fringsnfname: str,
) -> None:
    """Run FRING to calculate fring delays

    :param caldata: Calibrator data
    :type caldata: [type]
    :param snversion: TODO
    :type snversion: int
    :param clversion: TODO
    :type clversion: int
    :param sourcename: Name of source
    :type sourcename: str
    :param refant: Reference antenna
    :type refant: int
    :param fringsnfname: File name to save fring delays to
    :type fringsnfname: str
    """
    inttimesecs = 0.5  # Doesn't really matter if this is wrong
    applybandpasscal = False
    snrlimit = 6
    modeldata = None
    sumpols = False
    uvrange = [0, 0]
    zerorates = True
    vlbatasks.fring(
        caldata,
        snversion,
        clversion,
        SOLINTMINS,
        inttimesecs,
        sourcename,
        refant,
        applybandpasscal,
        snrlimit,
        SUMIFS,
        modeldata,
        sumpols,
        uvrange,
        zerorates,
        DELAYWINDOW,
        RATEWINDOW,
    )

    # Write SN table to disk
    if os.path.exists(fringsnfname):
        os.system("rm -f " + fringsnfname)
    vlbatasks.writetable(caldata, "SN", snversion, fringsnfname)


def correct_leakage(
    caldata,
    snversion: int,
    clversion: int,
    sourcename: str,
    refant: int,
    xpolmodelfile: str,
    xpolsnfname: str,
) -> None:
    """Correct for leakage if needed

    :param caldata: Calibrator data
    :type caldata: [type]
    :param snversion: TODO
    :type snversion: int
    :param clversion: TODO
    :type clversion: int
    :param sourcename: Source name
    :type sourcename: str
    :param refant: Reference antenna
    :type refant: int
    :param xpolmodelfile: X polarisation model file
    :type xpolmodelfile: str
    :param xpolsnfname: File containing X polarisation solutions
    :type xpolsnfname: str
    """
    xpolscan = 1
    if not os.path.exists(xpolmodelfile):
        print("Can't find xpol delay model  " + xpolmodelfile)
        print("Aborting!!")
        sys.exit(1)
    xpolmodel = AIPSImage("LKGSRC", "CLEAN", 1, 1)
    if xpolmodel.exists():
        xpolmodel.zap()
    vlbatasks.fitld_image(xpolmodelfile, xpolmodel)
    xpolsolintmins = 1
    inttimesecs = 0.5  # Doesn't matter if this is wrong
    if os.path.exists(xpolsnfname):
        os.remove(xpolsnfname)
    vlbatasks.xpoldelaycal(
        caldata,
        clversion,
        refant,
        sourcename,
        xpolscan,
        xpolmodel,
        xpolsolintmins,
        inttimesecs,
        xpolsnfname,
        DELAYWINDOW,
        RATEWINDOW,
    )
    vlbatasks.loadtable(caldata, xpolsnfname, snversion)
    vlbatasks.applysntable(caldata, snversion, "2PT", clversion, args.refant)


def run_bandpass(
    caldata,
    clversion: int,
    bpversion: int,
    sourcename: str,
    bpfname: str,
    cpasspoly: int,
    bpass: bool,
) -> None:
    """Run bandpass correction. Defaults to using CPASS unless --bpass
    is specified.

    :param caldata: Calibrator data
    :type caldata: [type]
    :param clversion: TODO
    :type clversion: int
    :param bpversion: TODO
    :type bpversion: int
    :param sourcename: Source name
    :type sourcename: str
    :param bpfname: Filename to save bandpass corrections to
    :type bpfname: str
    :param cpasspoly: Order of polynomial for CPASS
    :type cpasspoly: int
    :param bpass: If True, use BPASS instead of CPASS
    :type bpass: bool
    """
    scannumber = 1
    if bpass:
        vlbatasks.bpass(caldata, sourcename, clversion, scannumber)
    else:
        vlbatasks.cpass(
            caldata,
            sourcename,
            clversion,
            scannumber,
            None,
            cpasspoly,
        )

    # Write BP table to disk
    if os.path.exists(bpfname):
        os.system("rm -f " + bpfname)
    vlbatasks.writetable(caldata, "BP", bpversion, bpfname)


def plot_bandpass(caldata, bpversion: int, bptableplotfname: str) -> None:
    """Plot the bandpass corrections.

    :param caldata: Calibrator data
    :type caldata: [type]
    :param bpversion: TODO
    :type bpversion: int
    :param bptableplotfname: File to save plots into
    :type bptableplotfname: str
    """
    plotbptable = True
    vlbatasks.plotbandpass(
        caldata,
        bpversion,
        plotbptable,
        PLOTSPERPAGE,
        bptableplotfname,
    )


def run_selfcal(
    caldata,
    clversion: int,
    sourcename: str,
    refant: int,
    flux: float,
    selfcalsnfname: str,
) -> None:
    """Run selfcal

    :param caldata: Calibrator data
    :type caldata: [type]
    :param clversion: TODO
    :type clversion: int
    :param sourcename: Source name
    :type sourcename: str
    :param refant: Reference antenna
    :type refant: int
    :param flux: Calibrator flux in Jy
    :type flux: float
    :param selfcalsnfname: File to save selfcal solutions to
    :type selfcalsnfname: str
    """
    splitsnversion = 1
    doampcal = True
    dostokesi = False
    soltype = "L1R"
    selfcalsnr = 5
    splitcaldata = AIPSUVData(sourcename, OUTKLASS, 1, 1)
    if splitcaldata.exists():
        splitcaldata.zap()
    vlbatasks.split(caldata, clversion, OUTKLASS, sourcename)
    for i in range(1, 300):
        todeletedata = AIPSUVData(sourcename, "CALIB", 1, i)
        if todeletedata.exists():
            todeletedata.zap()
    vlbatasks.singlesource_calib(
        splitcaldata,
        flux,
        splitsnversion,
        refant,
        doampcal,
        SOLINTMINS,
        dostokesi,
        soltype,
        selfcalsnr,
        SUMIFS,
    )

    # Write SN table to disk
    if os.path.exists(selfcalsnfname):
        os.system("rm -f " + selfcalsnfname)
    vlbatasks.writetable(splitcaldata, "SN", 1, selfcalsnfname)


def plot_xcor(
    caldata,
    clversion: int,
    bpversion: int,
    xcorplotsmooth: int,
    uncalxcorplotfname: str,
    allcalxcorplotfname: str,
) -> None:
    """Plot calibrated and uncalibrated cross-correlation results

    :param caldata: Calibrator data
    :type caldata: [type]
    :param clversion: TODO
    :type clversion: int
    :param bpversion: TODO
    :type bpversion: int
    :param xcorplotsmooth: Length of the smoothing kernel in channels
        for xcor plotting
    :type xcorplotsmooth: int
    :param uncalxcorplotfname: File to save uncalibrated
        cross-correlation plots to
    :type uncalxcorplotfname: str
    :param allcalxcorplotfname: File to save calibrated
        cross-correlation plots to
    :type allcalxcorplotfname: str
    """
    plotbptable = False
    ifs = [0, 0]
    vlbatasks.plotbandpass(
        caldata,
        -1,
        plotbptable,
        PLOTSPERPAGE,
        uncalxcorplotfname,
        0,
        ifs,
        xcorplotsmooth,
    )
    vlbatasks.plotbandpass(
        caldata,
        bpversion,
        plotbptable,
        PLOTSPERPAGE,
        allcalxcorplotfname,
        clversion,
        ifs,
        xcorplotsmooth,
    )


def run_split(data, outfname: str, clversion: int, sourcename: str) -> None:
    """Run SPLIT on provided data

    :param data: Data to operate on
    :type data: [type]
    :param outfname: Filename to save output to
    :type outfname: str
    :param clversion: TODO
    :type clversion: int
    :param sourcename: Source name
    :type sourcename: str
    """
    outputdata = vlbatasks.zapAndCreateUVData(
        "CRAFTSRC", "SPLIT", AIPSDISK, SEQNO
    )
    vlbatasks.splitmulti(data, clversion, OUTKLASS, sourcename, SEQNO)
    vlbatasks.writedata(data, outfname + ".unaveraged", True)
    vlbatasks.writedata(outputdata, outfname, True)


def write_readme(
    calsolnfnames: "tuple[str]",
    calibrator: str,
    xpolmodelfile: str,
    reffreqs: "list[float]",
):
    """Write a README file for the calibration and a tarball with it and
    all the calibration solutions

    :param calsolnfnames: Tuple of calibration solution filenames as
        determined by `soln_fnames()`
    :type calsolnfnames: tuple[str]
    :param calibrator: Calibrator FITS file
    :type calibrator: str
    :param xpolmodelfile: File containing model used for xpol delay
        correction
    :type xpolmodelfile: str
    :param reffreqs: List of reference frequencies as determined by
        the calibrator data
    :type reffreqs: list[float]
    """
    bpfname = calsolnfnames[0]
    fringsnfname = calsolnfnames[1]
    selfcalsnfname = calsolnfnames[2]
    xpolsnfname = calsolnfnames[3]
    readmefname = calsolnfnames[7]
    calibtarballfile = calsolnfnames[8]

    readmeout = open(readmefname, "w")
    tarinputfiles = "{} {} {}".format(
        fringsnfname.split("/")[-1],
        selfcalsnfname.split("/")[-1],
        bpfname.split("/")[-1],
    )
    readmeout.write("This calibration was derived as follows:\n")
    readmeout.write("Calibrator file: %s\n" % calibrator)
    readmeout.write("Run on host: %s\n" % socket.gethostname())
    readmeout.write("At time: %s\n\n" % (str(datetime.datetime.now())))
    readmeout.write(
        "The following set of files was produced and used for calibration:\n"
    )
    readmeout.write(
        "%s (frequency-independent delay and phase from FRING)\n"
        % fringsnfname.split("/")[-1]
    )
    readmeout.write(
        "%s (frequency-independent complex gain [mostly just amplitude] from CALIB to set absolute flux scale)\n"
        % selfcalsnfname.split("/")[-1]
    )
    if xpolmodelfile != "":
        readmeout.write(
            "%s (frequency-independent, antenna-independent X-Y delay from FRING)\n"
            % xpolsnfname.split("/")[-1]
        )
        tarinputfiles = tarinputfiles + " " + xpolsnfname.split("/")[-1]
    readmeout.write(
        "%s (frequency-dependent complex gain from CPASS [polynomial bandpass fit])\n\n"
        % bpfname.split("/")[-1]
    )
    readmeout.write(
        "Remember that the delay specified in the SN tables generates zero phase at the reference frequency of the observation\n"
    )
    readmeout.write(
        'This reference frequency is set per subband (AIPS "IF"), but CRAFT datasets should only have one IF and hence one reference frequency.\n'
    )
    readmeout.write("Reference frequency(s) of this file:\n")
    for i, reffreq in enumerate(reffreqs):
        readmeout.write("AIPS IF %d ref (MHz): %.9f\n" % (i, reffreq))
    readmeout.write(
        "\nFinally I note that long-term, we really should also be solving for the leakage and writing both it and the parallactic angle corrections out.\n"
    )
    readmeout.close()
    if os.path.exists(calibtarballfile):
        os.system("rm -f " + calibtarballfile)
    os.system(f"tar cvzf {calibtarballfile} {readmefname} {tarinputfiles}")


if __name__ == "__main__":
    _main()
