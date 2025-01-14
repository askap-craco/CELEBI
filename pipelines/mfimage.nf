nextflow.enable.dsl=2

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
include {do_ref_correlation; do_correlation; get_start_mjd; difx_to_fits} from '../pipelines/correlate'
include {apply_offset} from '../pipelines/localise'

// Cards and FPGAs to be processed. Override these in a config file to cut out
// data. The lowest card-fpga pair is used as a reference correlation.
params.cards = ["1", "2", "3", "4", "5", "6", "7"]
cards = Channel.fromList(params.cards)
params.fpgas = ["0", "1", "2", "3", "4", "5"]
fpgas = Channel.fromList(params.fpgas)
card_fpgas = cards.combine(fpgas)
    .filter{ !(it[0] == params.cards.min() & it[1] == params.fpgas.min()) }
ref_card_fpga = cards.min().combine(fpgas.min())

// Other params
localise_dir = "$projectDir/../localise"
params.out_dir = "${params.publish_dir}/${params.label}"






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
        ml apptainer
        set -a
        set -o allexport

        # get args

        args="-i $data"
        args="\$args -s $summary"
        args="\$args -b $binconfig"
        args="\$args -p $polyco"
        args="\$args --thres $params.mf_thres"
        args="\$args --tN $params.mf_tN"
        args="\$args --rms_g $params.mf_rms_g"
        args="\$args --rms_w $params.mf_rms_w"
        args="\$args --rfi_w $params.mf_rfi_w"
        args="\$args --rfi_g $params.mf_rfi_g"

        if [ $params.mf_tw == 'true' ]; then
            args="\$args --tw"
        fi

        if [ $params.mf_fw == 'true' ]; then
            args="\$args --fw"
        fi


        # create new binconfig, polyco and wmask

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/make_mf_binconfig.py \$args'

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

    input:
        path corr_data
        path wmask

    output:
        path "D2DW/c*_f*", emit: mf_data
        path "*.png"

    script:
        """
        ml apptainer
        set -a
        set -o allexport

        # get args

        args="-d ./"
        args="\$args -w $wmask"
        args="\$args -p"


        # apply weights, rfi sub and scrunch correlated data

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/apply_mf_weights.py \$args'
        
        
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
                        binconfig.first(), polyco.first(), get_inttime.out.int_time.first(), 
                        startmjd, ref_correlation.combine(card_fpgas), fcm.first()).cx_fy

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

    input:
        path mf_fits
        path flux_cal_solns
    
    output:
        path "mf*"
        path "*calibrated_uv.ms"
        path "mf.jmfit", emit: mf_jmfit
    
    script:
        """
        # set up container stuff and tar flux cal solutions

        ml apptainer
        set -a
        set -o allexport
        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"
        
        tar -xzvf $flux_cal_solns


        # set up args
        args="--targetonly"
        args="\$args -t $mf_fits"
        args="\$args -r $params.refant"
        args="\$args -i"
        args="\$args -j"
        args="\$args --cleanmfs"
        args="\$args --pols=I"
        args="\$args --imagesize=$params.mfimagesize"
        args="\$args --pixelsize=$params.mfpixelsize"
        args="\$args --imagename=mf"
        args="\$args -a 16"
        aipsid="\$((RANDOM%8192))"
        args="\$args -u \$aipsid"
        args="\$args --src=$params.target"
        args="\$args --nmaxsources=1"
        args="\$args --findsourcescript=$localise_dir/get_pixels_from_field.py"
        args="\$args --findsourcescript2=$localise_dir/get_pixels_from_field2.py"

        if [ "$params.finderflagfile" != "" ] && [ "$params.finderflagfile" != "null" ]; then
            args="\$args --tarflagfile=$params.finderflagfile"
        fi


        # calibrate and image

        apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ParselTongue $localise_dir/calibrateFRB.py \$args' 

        rm -rf \$aips_dir  

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
        mf_fits = difx_to_fits("mf", mf_data.collect(), polyco, "mf").fits


        // Calibrate and image FRB
        mf_jmfit = mf_calibrate_and_image(mf_fits, flux_cal_solns).mf_jmfit


        // Apply RACS field source offset (if applicable)
        racs_offset = "${params.out_dir}/position/offset0.dat"
        racs_doffset = "${params.out_dir}/position/offsetfit.txt"

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

