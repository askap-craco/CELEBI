nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_pcal } from './correlate'
include { beamform as bform_pcal } from './beamform'
include { image_polcal; determine_pol_cal_solns as get_cal_pcal } from './calibration'
include { flag_proper as flagdat } from './flagging'

params.out_dir = "${params.publish_dir}/${params.label}"

workflow process_pol_cal {
    /*
        Process voltages to obtain polarisation calibration solutions

        Take
            flux_cal_solns: path
                Flux calibrator solutions tarball
        
        Emit
            pol_cal_solns: val/path
                Polarisation calibration solutions. Either an empty string if 
                solutions were not found (e.g. because the data has not been 
                flagged yet) or a text file containing the solutions.
    */    

    take:
        flux_cal_solns

    main:
        label = "${params.label}_polcal"
        empty_file = create_empty_file("file")

        // Correlation
        polcal_fits_path = "${params.out_dir}/loadfits/polcal/${params.label}_polcal.fits"
        if(new File(polcal_fits_path).exists()) {
            fits = Channel.fromPath(polcal_fits_path)
        }
        else {
            fits = corr_pcal(
                label, params.data_polcal, params.ra_polcal, params.dec_polcal, 
                empty_file, empty_file, empty_file, "polcal"
            ).fits
        }

        // Flagging
        if("${params.autoflag}" == "true") {
            polcal_fits_flagged = "${params.out_dir}/loadfits/polcal/${params.label}_polcal_f.fits"
        
            if(new File(polcal_fits_flagged).exists()) {
                    outfits = Channel.fromPath(polcal_fits_flagged)
            }
            else {
                outfits = flagdat(fits,polcal_fits_flagged, "cal").outfile
            }
            fits = outfits
        }

        // Calibration
        polcal_jmfit_path = "${params.out_dir}/polcal/polcal.jmfit"
        if(new File(polcal_jmfit_path).exists()) {
            pos = Channel.fromPath(polcal_jmfit_path)
        }
        else if(params.calibrate) {
                if(params.polflagfile == "") {
                    println "No polcal flag file!"
                    System.exit(1)
                }
                pos = image_polcal(
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
