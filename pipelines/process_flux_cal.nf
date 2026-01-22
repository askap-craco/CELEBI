nextflow.enable.dsl=2
nextflow.enable.strict=true // be less generous

include { create_empty_file } from './utils'
include { correlate as corr_fcal } from './correlate'
include { image_fluxcal; determine_flux_cal_solns as cal_fcal } from './calibration'
include { flag_proper as flagdat } from './flagging'


utils_dir    = "${projectDir}/../utils/"
beamform_dir = "${projectDir}/../beamform/"
localise_dir = "${projectDir}/../localise/"


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
            fits: val/path
                Flagged and calibrated fluxcal visibilities
    */
    take:
        fcm

    main:
        label = "${params.label}_fluxcal"
        
        fluxcal_fits_path = "${params.out_dir}/loadfits/fluxcal/${params.label}_fluxcal.fits"
        fluxcal_solns_path = "${params.out_dir}/fluxcal/calibration_noxpol_${params.target}.tar.gz"
        fcm_delayfix_path = "${params.out_dir}/fluxcal/fcm_delayfix.txt"
        
        // Correlation
        empty_file = create_empty_file("binconfig")
        
        if(params.binconfig_fluxcal == "") {
          binconfigpath = empty_file
          polycopath = empty_file
          inttimepath = empty_file
        }
        else {
          binconfigpath = file(params.binconfig_fluxcal)
          polycopath = file(params.polyco_fluxcal)
          inttimepath = file(params.inttime_fluxcal)
        }   
        
        fits = corr_fcal(
            label, params.data_fluxcal, params.ra_fluxcal, params.dec_fluxcal, 
            binconfigpath, polycopath, inttimepath, "fluxcal", fcm
        ).fits
        
		// Flagging
        if(!params.noflag) {
            fluxcal_fits_flagged = "${params.out_dir}/loadfits/fluxcal/${params.label}_fluxcal_f.fits"
            outfits = flagdat(fits,fluxcal_fits_flagged,"cal").outfile            
            fits = outfits
        }

		// Calibration         
        cal_fcal(fits, fcm)
        flux_cal_solns = cal_fcal.out.solns
        fcm_delayfix = cal_fcal.out.fcm_delayfix
        
    emit:
        flux_cal_solns
        fcm_delayfix
        fits
}
