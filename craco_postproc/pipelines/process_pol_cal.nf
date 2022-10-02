nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_pcal } from './correlate'
include { beamform as bform_pcal } from './beamform'
include { apply_flux_cal_solns_polcal as apply_cal_pcal; determine_pol_cal_solns as get_cal_pcal } from './calibration'
include { generate_binconfig } from './localise'

workflow process_pol_cal {
    /*
        Process voltages to obtain polarisation calibration solutions

        Take
            label: val
                FRB name and context as a string (no spaces)
            target: val
                FRB name as string (e.g. 181112)
            data: val
                Absolute path to pol cal data base directory (the dir. with the 
                ak* directories)
            data_frb: val
                Absolute path to FRB data base directory (the dir. with the ak*
                directories)
            snoopy: val
                Absolute path to snoopyv2.log file of FRB trigger
            fcm: val
                Absolute path to fcm (hardware delays) file
            ra0: val
                FRB right ascension initial guess as "hh:mm:ss"
            dec0: val
                FRB declination initial guess as "dd:mm:ss"
            polflagfile: val
                Absolute path to AIPS flag file for polarisation calibrator. If
                set to a blank string, the workflow will end before calibrating
            flux_cal_solns: path
                Flux calibrator solutions tarball
            cpasspoly: val
                Order of polynomial to fit bandpass with
            num_ints: val
                TODO: deprecate
            int_len: val
                TODO: deprecate
            offset: val
                TODO: deprecate
            dm: val
                DM to dedisperse to in pc/cm3
            centre_freq: val
                Central frequency of data in MHz
        
        Emit
            pol_cal_solns: val/path
                Polarisation calibration solutions. Either an empty string if 
                solutions were not found (e.g. because the data has not been 
                flagged yet) or a text file containing the solutions.
    */    

    take:
        label
        target
        data
        data_frb
        snoopy
        fcm
        ra0
        dec0
        polflagfile
        cpasspoly
        flux_cal_solns
        num_ints
        int_len
        offset
        dm
        centre_freq

    main:
        empty_file = create_empty_file("file")

        // Correlation
        polcal_fits_path = "${params.publish_dir}/${params.label}/loadfits/polcal/${params.label}_polcal.fits"
        if(new File(polcal_fits_path).exists()) {
            fits = Channel.fromPath(polcal_fits_path)
        }
        else {
            binconfig = generate_binconfig(data_frb, snoopy)
            fits = corr_pcal(
                label, data, fcm, ra0, dec0, empty_file, binconfig.polyco, 0, "polcal"
            )
        }

        // Calibration
        polcal_jmfit_path = "${params.publish_dir}/${params.label}/polcal/polcal.jmfit"
        if(new File(polcal_jmfit_path).exists()) {
            pos = Channel.fromPath(polcal_jmfit_path)
        }
        else if(params.calibrate) {
                if(params.polflagfile == "") {
                    println "No polcal flag file!"
                    System.exit(1)
                }
                pos = apply_cal_pcal(
                    fits, flux_cal_solns, polflagfile, target, cpasspoly
                ).jmfit
        }

        // Beamforming
        if(params.beamform) {
            bform_pcal(
                label, data, fcm, pos, flux_cal_solns, empty_file,
                num_ints, int_len, offset, dm, centre_freq, "-ds -IQUV"
            )
            pol_cal_solns = get_cal_pcal(bform_pcal.out.htr_data).pol_cal_solns
        }   
        else {
            pol_cal_solns = ""
        }

    emit:
        pol_cal_solns
}
