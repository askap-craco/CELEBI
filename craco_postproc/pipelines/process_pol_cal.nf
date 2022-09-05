nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_pcal } from './correlate'
include { beamform as bform_pcal } from './beamform'
include { apply_flux_cal_solns_polcal as apply_cal_pcal; determine_pol_cal_solns as get_cal_pcal } from './calibration'

workflow process_pol_cal {
    /*
        Process voltages to obtain polarisation calibration solutions

        Take
            polyco: path
                Output of generate_binconfig created from FRB data and snoopy
                log
            flux_cal_solns: path
                Flux calibrator solutions tarball
        
        Emit
            pol_cal_solns: val/path
                Polarisation calibration solutions. Either an empty string if 
                solutions were not found (e.g. because the data has not been 
                flagged yet) or a text file containing the solutions.
    */    

    take:
        polyco
        flux_cal_solns

    main:
        label = "${params.label}_polcal"
        empty_file = create_empty_file("file")

        // Correlation
        polcal_fits_path = "${params.publish_dir}/${params.label}/loadfits/polcal/${params.label}_polcal.fits"
        if(new File(polcal_fits_path).exists()) {
            fits = Channel.fromPath(polcal_fits_path)
        }
        else {
            fits = corr_pcal(
                label, params.data_polcal, params.ra_polcal, params.dec_polcal, 
                empty_file, polyco, 0, "polcal"
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
                    fits, flux_cal_solns, params.polflagfile
                ).jmfit
        }

        // Beamforming
        if(params.beamform) {
            bform_pcal(
                label, params.data_polcal, pos, flux_cal_solns, empty_file, 
                params.dm_polcal, params.centre_freq_polcal, "-ds -IQUV"
            )
            pol_cal_solns = get_cal_pcal(bform_pcal.out.htr_data).pol_cal_solns
        }   
        else {
            pol_cal_solns = ""
        }

    emit:
        pol_cal_solns
}
