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

# Global variables and option parsing
try:
    aipsver = os.environ["PSRVLBAIPSVER"]
except KeyError:
    aipsver = "31DEC18"


def _main():
    args = get_args()

    AIPS.userno = args.userno
    refant = args.refant
    xpolmodelfile = args.xpoldelaymodelfile
    xcorplotsmooth = args.xcorplotsmooth
    snversion = 1
    clversion = 1
    aipsdisk = 1

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

    # Load up the target data if needed
    if not args.calibrateonly:
        targetdata = vlbatasks.zapAndCreateUVData(
            "CRAFTTARG", "UVDATA", aipsdisk, 1
        )
        vlbatasks.fitld_corr(args.target, targetdata, [], "", 0.0001)
        if args.uvsrt:
            sortedtargetdata = vlbatasks.zapAndCreateUVData(
                "CRAFTTARG", "UVSRT", aipsdisk, 1
            )
            vlbatasks.uvsrt(targetdata, sortedtargetdata)
            targetdata.zap()
            targetdata = sortedtargetdata

        # Get the number of channels in the dataset
        numchannels = vlbatasks.getNumChannels(targetdata)

    # Load up the calibrator data
    if not args.targetonly:
        caldata = vlbatasks.zapAndCreateUVData(
            "CRAFTCAL", "UVDATA", aipsdisk, 1
        )
        vlbatasks.fitld_corr(args.calibrator, caldata, [], "", 0.0001)
        if args.uvsrt:
            sortedcaldata = vlbatasks.zapAndCreateUVData(
                "CRAFTCAL", "UVSRT", aipsdisk, 1
            )
            vlbatasks.uvsrt(caldata, sortedcaldata)
            caldata.zap()
            caldata = sortedcaldata

        # Get the reference frequency of the dataset
        reffreqs = []
        fqtable = caldata.table("FQ", 1)
        for row in fqtable:
            try:
                for iffreq in row.if_freq:
                    freqentry = float(iffreq) + float(caldata.header.crval[2])
                    reffreqs.append(
                        float(iffreq) + float(caldata.header.crval[2])
                    )
            except (AttributeError, TypeError):
                freqentry = float(row.if_freq) + float(caldata.header.crval[2])
                reffreqs.append(
                    float(row.if_freq) + float(caldata.header.crval[2])
                )

        # Flag the calibrator data, if desired
        if args.flagfile != "":
            vlbatasks.userflag(caldata, 1, args.flagfile)
        if args.shadow:
            shadowdiameter = float(args.shadow[0])
            xtalkbl = float(args.shadow[1])
            vlbatasks.shadowflag(caldata, 1, shadowdiameter, xtalkbl)
            print("Shadowing diameter: " + str(shadowdiameter))
            print("Cross-talk baseline: " + str(xtalkbl))

    # Flag the target data, if desired
    if not args.calibrateonly:
        if args.tarflagfile != "":
            vlbatasks.userflag(targetdata, 1, args.tarflagfile)
        if args.shadow:
            shadowdiameter = float(args.shadow[0])
            xtalkbl = float(args.shadow[1])
            vlbatasks.shadowflag(targetdata, 1, shadowdiameter, xtalkbl)
            print("Shadowing diameter: " + str(shadowdiameter))
            print("Cross-talk baseline: " + str(xtalkbl))

    # Run CLCOR to correct PANG if needed
    if xpolmodelfile != "":
        if not args.targetonly:
            vlbatasks.clcor_pang(caldata, clversion)
        if not args.calibrateonly:
            vlbatasks.clcor_pang(targetdata, clversion)
        clversion = clversion + 1

    # Run FRING
    if not args.targetonly:
        solintmins = 1  # Long enough that we just get one solutions
        inttimesecs = 0.5  # Doesn't really matter if this is wrong
        applybandpasscal = False
        snrlimit = 6
        sumifs = False
        modeldata = None
        sumpols = False
        uvrange = [0, 0]
        zerorates = True
        delaywindow = 0  # Search everything
        ratewindow = -1  # Don't search rates
        vlbatasks.fring(
            caldata,
            snversion,
            clversion,
            solintmins,
            inttimesecs,
            args.sourcename,
            refant,
            applybandpasscal,
            snrlimit,
            sumifs,
            modeldata,
            sumpols,
            uvrange,
            zerorates,
            delaywindow,
            ratewindow,
        )

        # Write SN table to disk
        if os.path.exists(fringsnfname):
            os.system("rm -f " + fringsnfname)
        vlbatasks.writetable(caldata, "SN", snversion, fringsnfname)

    # Load FRING SN table into the target
    if not args.calibrateonly:
        vlbatasks.loadtable(targetdata, fringsnfname, snversion)

    # Calibrate
    if not args.targetonly:
        vlbatasks.applysntable(caldata, snversion, "SELN", clversion, refant)
    if not args.calibrateonly:
        vlbatasks.applysntable(
            targetdata, snversion, "SELN", clversion, refant
        )
    snversion += 1
    clversion += 1

    # Correct for leakage if needed
    # leakagedopol = 0
    if xpolmodelfile != "":
        # First the xpoldelays
        if not args.targetonly:
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
                args.sourcename,
                xpolscan,
                xpolmodel,
                xpolsolintmins,
                inttimesecs,
                xpolsnfname,
                delaywindow,
                ratewindow,
            )
            vlbatasks.loadtable(caldata, xpolsnfname, snversion)
            vlbatasks.applysntable(
                caldata, snversion, "2PT", clversion, refant
            )
        if not args.calibrateonly:
            vlbatasks.loadtable(targetdata, xpolsnfname, snversion)
            vlbatasks.applysntable(
                targetdata, snversion, "2PT", clversion, refant
            )
        snversion += 1
        clversion += 1

        # Then the leakage
        # leakagefilename = os.getcwd() + "/leakage.an"
        # hasbptable = False
        # leakagemodel = xpolmodel
        # leakageacalmins = 1
        # leakagepcalmins = 1
        # leakagescan = 1
        # hasbptable = False
        # leakageoutputfile = os.getcwd() + '/' + options.sourcename + "_leakagecal_uv.fits"
        # leakageuvrange = [0,0]
        # leakageweightit = 0
        # vlbatasks.leakagecalc(caldata, options.sourcename, leakagemodel, leakagefilename,
        #            refant, leakageacalmins, leakagepcalmins, leakagescan, clversion,
        #            hasbptable, leakageoutputfile, leakageuvrange, leakageweightit)
        # vlbatasks.deletetable(caldata, "AN", 1)
        # vlbatasks.loadtable(caldata, leakagefilename, 1)
        # vlbatasks.deletetable(targetdata, "AN", 1)
        # vlbatasks.loadtable(targetdata, leakagefilename, 1)
        # leakagedopol = 2
        # print "Need to actually use leakagedopol below here - aborting!"
        # sys.exit()

    # Run BPASS
    # scannumber = 1
    # bpversion = 1
    # vlbatasks.bpass(caldata, options.sourcename, clversion, scannumber)

    # Run bandpass correction - default to CPASS, unless --bpass is specified
    scannumber = 1
    bpversion = 1
    if not args.targetonly:
        if args.bpass:
            vlbatasks.bpass(caldata, args.sourcename, clversion, scannumber)
        else:
            vlbatasks.cpass(
                caldata,
                args.sourcename,
                clversion,
                scannumber,
                None,
                args.cpasspoly,
            )

        # Write BP table to disk
        if os.path.exists(bpfname):
            os.system("rm -f " + bpfname)
        vlbatasks.writetable(caldata, "BP", bpversion, bpfname)

        # Plot the bandpass table
        if not args.skipplot:
            bptableplotfname = os.path.abspath(f"bptable{xpol_prefix}{src}.ps")
            plotsperpage = 4
            plotbptable = True
            vlbatasks.plotbandpass(
                caldata,
                bpversion,
                plotbptable,
                plotsperpage,
                bptableplotfname,
            )

    # Load up the bandpass to the target
    if not args.calibrateonly:
        vlbatasks.loadtable(targetdata, bpfname, bpversion)

    # Run selfcal
    outklass = "SPLIT"
    if not args.targetonly:
        applybandpasscal = True
        splitsnversion = 1
        doampcal = True
        dostokesi = False
        soltype = "L1R"
        selfcalsnr = 5
        splitcaldata = AIPSUVData(args.sourcename, outklass, 1, 1)
        if splitcaldata.exists():
            splitcaldata.zap()
        vlbatasks.split(caldata, clversion, outklass, args.sourcename)
        for i in range(1, 300):
            todeletedata = AIPSUVData(args.sourcename, "CALIB", 1, i)
            if todeletedata.exists():
                todeletedata.zap()
        vlbatasks.singlesource_calib(
            splitcaldata,
            args.flux,
            splitsnversion,
            args.refant,
            doampcal,
            solintmins,
            dostokesi,
            soltype,
            selfcalsnr,
            sumifs,
        )

        # Write SN table to disk
        if os.path.exists(selfcalsnfname):
            os.system("rm -f " + selfcalsnfname)
        vlbatasks.writetable(splitcaldata, "SN", 1, selfcalsnfname)

    # Load up the selfcal SN table
    if not args.calibrateonly:
        vlbatasks.loadtable(targetdata, selfcalsnfname, snversion)
    if not args.targetonly:
        vlbatasks.loadtable(caldata, selfcalsnfname, snversion)

    # Calibrate
    if not args.targetonly:
        vlbatasks.applysntable(caldata, snversion, "SELN", clversion, refant)
    if not args.calibrateonly:
        vlbatasks.applysntable(
            targetdata, snversion, "SELN", clversion, refant
        )
    snversion += 1
    clversion += 1

    # Plot the uncalibrated and calibrated cross-correlation results if desired
    if not args.targetonly:
        if not args.skipplot:
            uncalxcorplotfname = os.path.abspath(
                f"uncalxcor{xpol_prefix}{src}.ps"
            )
            allcalxcorplotfname = os.path.abspath(
                f"allcalxcor{xpol_prefix}{src}.ps"
            )
            plotbptable = False
            plotsperpage = 4
            ifs = [0, 0]
            vlbatasks.plotbandpass(
                caldata,
                -1,
                plotbptable,
                plotsperpage,
                uncalxcorplotfname,
                0,
                ifs,
                xcorplotsmooth,
            )
            vlbatasks.plotbandpass(
                caldata,
                bpversion,
                plotbptable,
                plotsperpage,
                allcalxcorplotfname,
                clversion,
                ifs,
                xcorplotsmooth,
            )

    # Run SPLIT and write output data for calibrator
    seqno = 1
    if not args.targetonly:
        outputdata = vlbatasks.zapAndCreateUVData(
            "CRAFTSRC", "SPLIT", aipsdisk, seqno
        )
        vlbatasks.splitmulti(
            caldata, clversion, outklass, args.sourcename, seqno
        )
        vlbatasks.writedata(caldata, caloutfname + ".unaveraged", True)
        vlbatasks.writedata(
            outputdata, caloutfname, True
        )  # TODO: Make this optional

    # Run SPLIT and write output data for target
    seqno = 1
    if not args.calibrateonly:
        outputdata = vlbatasks.zapAndCreateUVData(
            "CRAFTSRC", "SPLIT", aipsdisk, seqno
        )
        vlbatasks.splitmulti(
            targetdata, clversion, outklass, args.sourcename, seqno
        )
        vlbatasks.writedata(
            targetdata, targetoutfname + ".unaveraged", True
        )  # TODO: Make this optional
        vlbatasks.writedata(outputdata, targetoutfname, True)

    # Create a README file for the calibration and a tarball with it plus all the calibration
    if not args.targetonly:
        readmeout = open(f"README{xpol_prefix}{src}.calibration", "w")
        tarinputfiles = "{} {} {}".format(
            fringsnfname.split("/")[-1],
            selfcalsnfname.split("/")[-1],
            bpfname.split("/")[-1],
        )
        readmeout.write("This calibration was derived as follows:\n")
        readmeout.write("Calibrator file: %s\n" % args.calibrator)
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
        calibtarballfile = f"calibration{xpol_prefix}{src}.tar.gz"
        if os.path.exists(calibtarballfile):
            os.system("rm -f " + calibtarballfile)
        os.system(
            f"tar cvzf {calibtarballfile} README{xpol_prefix}{src}.calibration {tarinputfiles}"
        )

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


def do_imaging(args: argparse.Namespace) -> None:
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
                os.system("casa --nologger -c imagescript.py")
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


def soln_fnames(
    xpolmodelfile: str, src: str
) -> "tuple[str, str, str, str, str, str, str]":
    """Determine filenames for the calibration solutions.

    :param xpolmodelfile: X polarisation delay model filename
    :type xpolmodelfile: str
    :param src: Source name (i.e. FRB/Vela/etc.)
    :type src: str
    :return: Filenames for:
        - Bandpasses
        - Fring delays
        - Selcal solutions
        - X polarisation fring delays
        - Bandpass table plots
        - Uncalibrated cross-correlation plots
        - Calibrated cross-correlation plots
    :rtype: tuple[str, str, str, str, str, str, str]
    """
    if xpolmodelfile != "":
        xpol_prefix = "_xpol"
    else:
        xpol_prefix = "_noxpol"
    if src != "":
        src = "_" + src
    else:
        src = src
    bpfname = os.path.abspath(f"bandpasses{xpol_prefix}{src}.bp.txt")
    fringsnfname = os.path.abspath(f"delays{xpol_prefix}{src}.sn.txt")
    selfcalsnfname = os.path.abspath(f"selfcal{xpol_prefix}{src}.sn.txt")
    xpolsnfname = os.path.abspath(f"xpolfring{xpol_prefix}{src}.sn.txt")
    bptableplotfname = os.path.abspath(f"bptable{xpol_prefix}{src}.ps")
    uncalxcorplotfname = os.path.abspath(f"uncalxcor{xpol_prefix}{src}.ps")
    allcalxcorplotfname = os.path.abspath(f"allcalxcor{xpol_prefix}{src}.ps")

    return (
        bpfname,
        fringsnfname,
        selfcalsnfname,
        xpolsnfname,
        bptableplotfname,
        uncalxcorplotfname,
        allcalxcorplotfname,
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


if __name__ == "__main__":
    _main()
