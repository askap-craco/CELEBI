nextflow.enable.dsl=2
nextflow.enable.strict=true // be less generous

include { create_empty_file } from './utils'
include { correlate as corr_pcal } from './correlate'
include { beamform as bform_pcal; gen_dspec as pcal_dspec } from './beamform'
include { image_polcal; determine_pol_cal_solns as get_cal_pcal } from './calibration'
include { flag_proper as flagdat } from './flagging'


utils_dir    = "${projectDir}/../utils/"
beamform_dir = "${projectDir}/../beamform/"
localise_dir = "${projectDir}/../localise/"


workflow process_pol_cal {
    /*
        Process voltages to obtain polarisation calibration solutions

        Take
            flux_cal_solns: path
                Flux calibrator solutions tarball
            fcm: path
                fcm file to use in correlation
            antspcal: val
                List of antennas
        
        Emit
            pol_cal_solns: val/path
                Polarisation calibration solutions. Either an empty string if 
                solutions were not found (e.g. because the data has not been 
                flagged yet) or a text file containing the solutions.
    */    

    take:
        flux_cal_solns
        fcm
        antspcal

    main:
        label = "${params.label}_polcal"
        empty_file = create_empty_file("file")
		
		polcal_fits_path = "${params.out_dir}/loadfits/polcal/${params.label}_polcal.fits" 
        if( params.localize || params.corrpcal ) {
            // Correlation          
            fits = corr_pcal(
                label, params.data_polcal, params.ra_polcal, params.dec_polcal, 
                empty_file, empty_file, empty_file, "polcal", fcm
            ).fits
        }
        else {
            fits = file(polcal_fits_path)
        }
        
        polcal_jmfit_path = "${params.out_dir}/polcal/polcal.jmfit" 
        if( params.localize || params.imgpcal ) {   
            // Flagging
            if(!params.noflag) {
                polcal_fits_flagged = "${params.out_dir}/loadfits/polcal/${params.label}_polcal_f.fits"            
                outfits = flagdat(fits,polcal_fits_flagged, "cal").outfile                
                fits = outfits
            }
            
            pos = image_polcal(
                fits, flux_cal_solns, params.polflagfile
            ).jmfit            
        }
        else {
            pos = file(polcal_jmfit_path)
        }

		// Beamforming
        polcal_solns_path = "${params.out_dir}/polcal/${params.label}_polcal_solutions.txt"
        xypath = "${params.out_dir}/htr/${label}_*_t_${params.dm_polcal}.npy"
        if( params.gethtr || params.beampcal ) {
            bform_pcal(
                label, params.data_polcal, pos, flux_cal_solns, empty_file, 
                params.dm_polcal, params.centre_freq_polcal,
                antspcal, fcm, "NONE"
            )
            xy = bform_pcal.out.xy
        }
        else {
            xy = file(xypath)
        }
        
        // Calculate polcal solutions
        if( params.gethtr || params.calcpcal ) {
            pcal_dspec(
                label, xy, params.dm_polcal, params.centre_freq_polcal
            )
            pol_cal_solns = get_cal_pcal(pcal_dspec.out.htr_data).pol_cal_solns
        }
        else {
            pol_cal_solns = file(polcal_solns_path)
        }
        
    emit:
        pol_cal_solns
}





