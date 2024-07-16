nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_fcal } from './correlate'
include { determine_flux_cal_solns as cal_fcal } from './calibration'
include { flag_proper as flagdat } from './flagging'

params.out_dir = "${params.publish_dir}/${params.label}"

workflow process_flux_cal {
    /*
        Process voltages to obtain flux+phase calibration solutions
        
        Take
            fcm: path
                fcm file to update. If an empty file, won't be updated

        Emit
            flux_cal_solns: val/path
                Flux calibration solutions. Either an empty string if solutions
                were not found (e.g. because the data has not been flagged yet)   
                or a tarball containing the solutions.
            fcm_delayfix: val/path
                Delayfixed FCM file. Either an empty string if not done 
                (e.g. because the data has not been flagged yet) or a txt file
    */
    take:
        fcm

    main:
        label = "${params.label}_fluxcal"
        // Correlation
        fluxcal_fits_path = "${params.out_dir}/loadfits/fluxcal/${params.label}_fluxcal.fits"
        // if(new File(fluxcal_fits_path).exists()) {
        //     fits = Channel.fromPath(fluxcal_fits_path)
        // }
        // else {
        empty_file = create_empty_file("binconfig")
        if(params.binconfig_fluxcal == "") {
          binconfigpath = empty_file
          polycopath = empty_file
          inttimepath = empty_file
        }
        else {
          binconfigpath = Channel.fromPath(params.binconfig_fluxcal) //.first()
          polycopath = Channel.fromPath(params.polyco_fluxcal) //.first()
          inttimepath = Channel.fromPath(params.inttime_fluxcal) //.first()
        }
        fits = corr_fcal(
            label, params.data_fluxcal, params.ra_fluxcal, params.dec_fluxcal, 
            binconfigpath, polycopath, inttimepath, "fluxcal", fcm
        ).fits
        // }
        fits.view()
        // Flagging
        if(!params.noflag) {
            fluxcal_fits_flagged = "${params.out_dir}/loadfits/fluxcal/${params.label}_fluxcal_f.fits"
        
            // if(new File(fluxcal_fits_flagged).exists()) {
            //         outfits = Channel.fromPath(fluxcal_fits_flagged)
            // }
            // else {
            outfits = flagdat(fits,fluxcal_fits_flagged,"cal").outfile
            // }
            
            fits = outfits
        }

        // Calibration
        fluxcal_solns_path = "${params.out_dir}/fluxcal/calibration_noxpol_${params.target}.tar.gz"
        fcm_delayfix_path = "${params.out_dir}/fluxcal/fcm_delayfix.txt"
        // if(new File(fluxcal_solns_path).exists() and new File(fcm_delayfix_path).exists()) {
        //     flux_cal_solns = Channel.fromPath(fluxcal_solns_path)
        //     fcm_delayfix = Channel.fromPath(fcm_delayfix_path)
        // }
        // else if(params.calibrate) {
        if(params.calibrate) {
            cal_fcal(fits, params.fluxflagfile, fcm)
            flux_cal_solns = cal_fcal.out.solns
            fcm_delayfix = cal_fcal.out.fcm_delayfix
        }
        else {
            flux_cal_solns = ""
            fcm_delayfix = ""
        }
        
    emit:
        flux_cal_solns
        fcm_delayfix
}
