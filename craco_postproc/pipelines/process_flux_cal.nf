nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_fcal } from './correlate'
include { determine_flux_cal_solns as cal_fcal } from './calibration'
 
workflow process_flux_cal {
    /*
        Process voltages to obtain flux+phase calibration solutions

        Take
            label: val
                FRB name and context as a string (no spaces)
            target: val
                FRB name as string (e.g. 181112)
            data: val
                Absolute path to flux calibrator data base directory (the dir. 
                with the ak* directories)
            binconfig: paths
                Output of generate_binconfig created from FRB data and snoopy
                log
            fcm: val
                Absolute path to fcm (hardware delays) file
            ra: val
                Flux calibrator right ascension as "hh:mm:ss"
            dec: val
                Flux calibrator declination as "dd:mm:ss"
            fluxflagfile: val
                Absolute path to AIPS flag file for flux calibrator. If set to
                a blank string, the workflow will end before calibrating.
            cpasspoly: val
                Order of polynomial to fit bandpass with
        
        Emit
            flux_cal_solns: val/path
                Flux calibration solutions. Either an empty string if solutions
                were not found (e.g. because the data has not been flagged yet)   
                or a tarball containing the solutions.
    */
    take:
        label
        target
        data
        binconfig
        fcm
        ra
        dec
        fluxflagfile
        cpasspoly

    main:
        // Correlation
        fluxcal_fits_path = "${params.publish_dir}/${params.label}/loadfits/fluxcal/${params.label}_fluxcal.fits"
        if(new File(fluxcal_fits_path).exists()) {
            fits = Channel.fromPath(fluxcal_fits_path)
        }
        else {
            empty_binconfig = create_empty_file("binconfig")
            fits = corr_fcal(
                label, data, fcm, ra, dec, empty_binconfig, binconfig.polyco, 0, 
                "fluxcal"
            )
        }

        // Calibration
        fluxcal_solns_path = "${params.publish_dir}/${params.label}/fluxcal/calibration_noxpol_${target}.tar.gz"
        if(new File(fluxcal_solns_path).exists()) {
            flux_cal_solns = Channel.fromPath(fluxcal_solns_path)
        }
        else if(params.calibrate) {
            if(params.fluxflagfile == "") {
                println "No fluxcal flag file!"
                System.exit(1)
            }
            flux_cal_solns = cal_fcal(fits, fluxflagfile, target, cpasspoly).solns
        }
        else {
            flux_cal_solns = ""
        }
        
    emit:
        flux_cal_solns
}
