nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_pcal } from './correlate'
include { beamform as bform_pcal } from './beamform'
include { image_polcal; determine_pol_cal_solns as get_cal_pcal } from './calibration'
include { flag_proper as flagdat } from './flagging'

params.out_dir = "${params.publish_dir}/${params.label}"
params.nants_pcal = params.nants

workflow process_pol_cal {
    /*
        Process voltages to obtain polarisation calibration solutions

        Take
            flux_cal_solns: path
                Flux calibrator solutions tarball
            fcm: path
                fcm file to use in correlation
        
        Emit
            pol_cal_solns: val/path
                Polarisation calibration solutions. Either an empty string if 
                solutions were not found (e.g. because the data has not been 
                flagged yet) or a text file containing the solutions.
    */    

    take:
        flux_cal_solns
        fcm

    main:
        label = "${params.label}_polcal"
        
        empty_file = create_empty_file("file")
        polcal_fits_path = "${params.out_dir}/loadfits/polcal/${params.label}_polcal.fits" 
        if( params.makeimage || params.corrcal ) {
            // Correlation          
            fits = corr_pcal(
                label, params.data_polcal, params.ra_polcal, params.dec_polcal, 
                empty_file, empty_file, empty_file, "polcal", fcm
            ).fits
        }
        else {
            fits = Channel.fromPath(polcal_fits_path)
        }
        
        polcal_jmfit_path = "${params.out_dir}/polcal/polcal.jmfit" 
        if( params.makeimage || params.impcal ) {   
            // Flagging
            if(!params.noflag) {
                polcal_fits_flagged = "${params.out_dir}/loadfits/polcal/${params.label}_polcal_f.fits"            
                outfits = flagdat(fits,polcal_fits_flagged, "cal").outfile                
                fits = outfits
            }

            // Calibration
            pos = image_polcal(
                fits, flux_cal_solns, params.polflagfile
            ).jmfit            
        }
        else {
            pos = Channel.fromPath(polcal_jmfit_path)
        }

        // Beamforming
        polcal_solns_path = "${params.out_dir}/polcal/${params.label}_polcal_solutions.txt"
        if( params.beamform || params.beampcal ) {
            bform_pcal(
                label, params.data_polcal, pos, flux_cal_solns, empty_file, 
                params.dm_polcal, params.centre_freq_polcal,
                params.nants_pcal, fcm, "NONE"
            )
            pol_cal_solns = get_cal_pcal(bform_pcal.out.htr_data).pol_cal_solns
        }   
        else {
            pol_cal_solns = Channel.fromPath(polcal_solns_path)
        }
        
    emit:
        pol_cal_solns
}
