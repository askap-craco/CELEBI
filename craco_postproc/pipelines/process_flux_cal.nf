nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_fcal } from './correlate'
include { determine_flux_cal_solns as cal_fcal } from './calibration'

workflow process_flux_cal {
    /*
        Process voltages to obtain flux+phase calibration solutions
        
        Emit
            flux_cal_solns: val/path
                Flux calibration solutions. Either an empty string if solutions
                were not found (e.g. because the data has not been flagged yet)   
                or a tarball containing the solutions.
    */
    main:
        label = "${params.label}_fluxcal"

        // Correlation
        fluxcal_fits_path = "${params.publish_dir}/${params.label}/loadfits/fluxcal/${params.label}_fluxcal.fits"
        if(new File(fluxcal_fits_path).exists()) {
            fits = Channel.fromPath(fluxcal_fits_path)
        }
        else {
            empty_file = create_empty_file("binconfig")
            fits = corr_fcal(
                label, params.data_fluxcal, params.ra_fluxcal, params.dec_fluxcal, 
                empty_file, empty_file, 0, "fluxcal"
            ).fits
        }

        // Calibration
        fluxcal_solns_path = "${params.publish_dir}/${params.label}/fluxcal/calibration_noxpol_${params.target}.tar.gz"
        if(new File(fluxcal_solns_path).exists()) {
            flux_cal_solns = Channel.fromPath(fluxcal_solns_path)
        }
        else if(params.calibrate) {
            if(params.fluxflagfile == "") {
                println "No fluxcal flag file!"
                System.exit(1)
            }
            flux_cal_solns = cal_fcal(fits, params.fluxflagfile).solns
        }
        else {
            flux_cal_solns = ""
        }
        
    emit:
        flux_cal_solns
}
