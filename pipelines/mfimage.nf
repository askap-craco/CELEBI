nextflow.enable.dsl=2   // Enable DSL2
nextflow.enable.strict=true // be less generous

/*
Author: Tyson Dial
Email: tdial@swin.edu.au
Created: 01/10/2024
Last Updated: 01/10/2024

Info:
    processes and workflow to apply a matched filter binconfig file to the
    correlated data when imaging.

*/


// Include nessesary scripts for correlating data (and other scripts)
include {do_ref_correlation; do_correlation; get_start_mjd; difx_to_fits} from './correlate'
include {apply_offset} from './localise'

flagging_dir = "${projectDir}/../flagging/"
utils_dir    = "${projectDir}/../utils/"
beamform_dir = "${projectDir}/../beamform/"
localise_dir = "${projectDir}/../localise/"


cards = Channel.fromList(params.cards)
fpgas = Channel.fromList(params.fpgas)
card_fpgas = cards.combine(fpgas)
    .filter{ !(it[0] == params.cards.min() & it[1] == params.fpgas.min()) }
ref_card_fpga = cards.min().combine(fpgas.min())



process generate_mf_binconfig {
    /*

        Create matched filter binconfig and polyco files for imaging frb.
        Also make wmask structure with time(freq-)dependent weights to apply to correlated
        data for improved imaging.

        Input:
            data: path
                Path to Stokes I HTR dynamic spectrum
            binconfig: path
                Path to binconfig file
            polyco: path
                Path to polyco file
            summary: path
                Path to summary file
            
        Output:
            binconfig: path
                New binconfig file
            polyco: path
                New polyco file
            wmask: path
                wmask struct holding time(freq-)dependent weights to apply to correlated 
                data
            *.png: path
                diagnostic plots 
            *.npy: path
                data files of diagnostic plots
            *.txt: path
                summary info paste-bin when creating new mf binconfig file
        

    */
    publishDir "${params.publish_dir}/${params.label}/mf", mode: "copy"
	
	label 'celebi'
	
    input:
        path data
        path binconfig
        path polyco
        path summary
    
    output:
        path "mf.binconfig", emit: binconfig
        path "mf.polyco", emit: polyco
        path "mf_wmask.npy", emit: wmask
        path "*.png"
        path "*.npy"
        path "*.txt"

    script:
        """
        source /opt/setup_proc_container 
        set -xu
		
        args=''
        if [ $params.mf_tw == 'true' ]; then
            args="\$args --tw"
        fi

        if [ $params.mf_fw == 'true' ]; then
            args="\$args --fw"
        fi
		# create new binconfig, polyco and wmask
		
        python3 $localise_dir/make_mf_binconfig.py \
        		-i $data \
        		-s $summary \
        		-b $binconfig \
        		-p $polyco \
        		--thres $params.mf_thres \
        		--tN $params.mf_tN \
        		--rms_g $params.mf_rms_g \
        		--rms_w $params.mf_rms_w \
        		--rfi_w $params.mf_rfi_w \
        		--rfi_g $params.mf_rfi_g \
                \$args

        """


}






process apply_mf_weights {
    /*
        Apply matched filter weights to correlated data and scrunch to a single bin

        Input:
            corr_data: path
                c*_f* directories containing .difx data (correlated data). This is only here
                to get links to the data
            wmask: path
                path to wmask structure containing time(freq-)dependent weights and finder/rfi
                masks
        
        Output:
            mf_data: path
                new c*_f* directories containing a single .difx bin (.b0000) which is
                the scrunched and weighted correlated data
            *.png: path
                diagnostic plots
            
    */
    publishDir "${params.publish_dir}/${params.label}/mf/weighted", mode: "copy"
	
	label 'celebi'
	
    input:
        path corr_data
        path wmask

    output:
        path "D2DW/c*_f*", emit: mf_data
        path "*.png"

    script:
        """
        source /opt/setup_proc_container 
        set -xu

       	# apply weights, rfi sub and scrunch correlated data
       	
       	python3 $localise_dir/apply_mf_weights.py \
       			-d ./ \
       			-w $wmask \
       			-p        
        
        """


}








process get_inttime {

    input:
        val int_time

    output:
        path "inttime.txt", emit: int_time
    
    script:
        """
        echo $int_time >> inttime.txt
        
        """
}






workflow correlate_frb {
    /*
        Correlate data using a binconfig and polyco file

        Take
            binconfig: path
                path to binconfig file
            polyco: path
                path to polyco file
            fcm: path
                fcm to use, ideally with delayfix 
    */
    take:
        binconfig
        polyco
        fcm
    main:

        // get int_time
        get_inttime(1.3824)

        // Get start mjd
        startmjd = get_start_mjd(params.data_frb)

        // Reference correlation
        ref_correlation = do_ref_correlation(params.label, params.data_frb, params.ra_frb, params.dec_frb, 
                        binconfig, polyco, get_inttime.out.int_time, startmjd, ref_card_fpga,
                        fcm).cx_fy


        // Do rest of correlations
        correlated_data = do_correlation(params.label, params.data_frb, params.ra_frb, params.dec_frb, 
                        binconfig, polyco, get_inttime.out.int_time, 
                        startmjd, ref_correlation.combine(card_fpgas), fcm).cx_fy

        // Combine correlations
        all_correlations = ref_correlation.concat(correlated_data).collect()
    
    emit:
        corr_data = all_correlations


}





process mf_calibrate_and_image {
    /*
        Calibrate and image .FITS mf finder bin

        Input:
            mf_fits: path
                .FITS file to image
            flux_cal_solns: path
                flux calibration solutions


    */

    publishDir "${params.publish_dir}/${params.label}/mf/image", mode: "copy"
	
	label 'celebi'
	label 'aips_tempfs'
	
    input:
        path mf_fits
        path flux_cal_solns
    
    output:
        path "mf*"
        path "*calibrated_uv.ms.tar"
        path "mf.jmfit", emit: mf_jmfit
    
    script:
        """
        # set up container stuff and tar flux cal solutions

        source /opt/setup_proc_container 
        set -xu
        
        export PATH=\$PATH:$params.casapath

        aipsid="\$((RANDOM%8192))"

        cp $mf_fits /JOBFS/.
        cp $flux_cal_solns /JOBFS/.
        cd /JOBFS
        
        tar -xzvf $flux_cal_solns
        
        mf_fits=$mf_fits
		
		if [ "$params.finderflagfile" != "" ] && [ "$params.finderflagfile" != "null" ]; then
            args=" --tarflagfile=$params.finderflagfile"
        else
            args=""
        fi
		
        export LC_CTYPE=C
        export LC_ALL=C
        export LANGUAGE=C
        
		ParselTongue $localise_dir/calibrateFRB.py \
				--targetonly \
				-t \$mf_fits \
                --maskpeakonly \
				-r $params.refant \
				-i \
				-j \
				--cleanmfs \
				--pols=I \
				--imagesize=$params.mfimagesize \
				--pixelsize=$params.mfpixelsize \
				--imagename=mf \
				-a 16 \
				-u \$aipsid \
				--src=$params.target \
				--nmaxsources=1 \
				--findsourcescript=$localise_dir/get_pixels_from_field.py \
				--findsourcescript2=$localise_dir/get_pixels_from_field2.py \
				\$args
		
		tar -cvf ${mf_fits}_calibrated_uv.ms.tar \${mf_fits%.fits}_calibrated_uv.ms
        cd - 

        cp -r /JOBFS/${mf_fits}_calibrated_uv.ms.tar .
		
        cp -r /JOBFS/mf.image .
        cp /JOBFS/mf.jmfit .		
		
        """


}




process if_racs {

    input:
        path offset_file
    output:
        env racs_exists, emit: racs_exists
    script:
        """
        racs_exists="0"
        if [ -f $offset_file ]; then
            racs_exists="1"
        fi
        
        """

}





workflow mf_image {
    /*
        Workflow to image FRB using a matched filter binconfig file with time(freq-)dependent
        weighting and rfi subtraction.

        Take
            stk_i: path
                path to Stokes I HTR dynamic spectra .npy file
            binconfig: path 
                path to original binconfig file
            polyco: path
                path to original polyco file
            summary: path
                path to summary file created at the end of an FRB run
            flux_cal_solns: 
                flux calibration solutions
            fcm:
                fcm to use, ideally with delayfix

    */
    take:
        stk_i
        binconfig
        polyco
        summary
        flux_cal_solns
        fcm
    main:
        // MAIN SCRIPT

        // Create matched filter binconfig, polyco and wmask files for time(freq-)dependent
        // weighting.
        mf_files = generate_mf_binconfig(stk_i, binconfig, polyco, summary)


        // Correlate FRB data using new binconfig and polyco files
        corr_data = correlate_frb(mf_files.binconfig, mf_files.polyco, fcm)


        // apply time(freq-)dependent weights, rfi subtraction and scrunch to a single
        // finder bin
        mf_data = apply_mf_weights(corr_data, mf_files.wmask).mf_data


        // convert .difx file to .FITS file
        mf_fits = difx_to_fits("${params.label}_mf", mf_data.collect(), polyco, "mf").fits


        // Calibrate and image FRB
        mf_jmfit = mf_calibrate_and_image(mf_fits, flux_cal_solns).mf_jmfit


        // Apply RACS field source offset (if applicable)
        racs_offset = file("${params.out_dir}/position/offset0.dat")
        racs_doffset = file("${params.out_dir}/position/offsetfit.txt")

        // first check "offset0.dat" file exists
        if_racs(racs_offset)

        if (if_racs.out.racs_exists.toInteger()) {

            // apply racs
            final_res = apply_offset(racs_offset, racs_doffset, mf_jmfit)
            mf_final_position = final_res.final_position

        }
        else {
            mf_final_position = Channel.empty()
        }
    
    emit:
        mf_final_position = mf_final_position

}

