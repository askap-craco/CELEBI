/*
    Processes and workflow related to beamforming and producing high-time 
    resolution data products
*/

nextflow.enable.dsl=2
nextflow.enable.strict=true // be less generous

include { get_start_mjd } from './correlate'
include { apply_pol_cal_solns } from './calibration'
include { filter_antenna } from './utils'


utils_dir    = "${projectDir}/../utils/"
beamform_dir = "${projectDir}/../beamform/"
localise_dir = "${projectDir}/../localise/"


polarisations = Channel
    .fromList(params.pols)


process create_calcfiles {
    /*
        Create the files that contain the delays required for 
        beamforming on a particular position.

        Input
            label: val
                FRB name and context of process instance as a string (no 
                spaces)
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
            pos: path
                File containing JMFIT statistics of position to beamform on
            fcm: path
                fcm to use
        
        Output
            Delay files: tuple(path, path)
                The .im and .calc files that are used by craftcor_tab.py to 
                calculate geometric delays for beamforming
    */

    label 'celebi'

    input:
        val label
        val data
        path pos
        path fcm

    output:
        tuple path("c1_f0/craftfrb*.im"), path("c1_f0/craftfrb*.calc")

    script:
        """
        source /opt/setup_proc_container
        set -xu

        startmjd=`python3 $localise_dir/get_start_mjd.py $data` 

        export CRAFTCATDIR="."

        clabel=$label
        contxt=\${clabel:(-6)}
        echo "\$contxt"

        if [ "$params.usepos" = "true" ] && [ "\$contxt" != "polcal" ]; then
            ra=$params.ra_frb
            dec=$params.dec_frb
            echo "Using given position \$ra \$dec"
        else
            ra=\$(grep "Actual RA" $pos)
            ra=\${ra:22}
            dec=\$(grep "Actual Dec" $pos)
            dec=\${dec:22}
            echo "Using JMFIT position \$ra \$dec"
        fi

        # Run processTimeStep.py with the --calconly flag to stop once calcfile is 
        # written
        # This is a clunky if statement to handle TLEs if needed, I'm sure this could be streamlined
        if [ "$params.tle_file" != "" ] && [ "$params.tle_file" != "null" ]; then
          python3 $localise_dir/processTimeStep.py -t=$data \
              --ra=\$ra \
              --dec=\$dec \
              --tlefile=$params.tle_file \
              --tleobject=$params.tle_object \
              -f=$fcm \
              -b=$params.nbits \
              --card=1 \
              -k \
              --name=$label \
              -o . \
              --freqlabel c1_f0 \
              --dir=$projectDir/../difx \
              --calconly \
              --startmjd=\$startmjd
        else
          python3 $localise_dir/processTimeStep.py -t=$data \
              --ra=\$ra \
              --dec=\$dec \
              -f=$fcm \
              -b=$params.nbits \
              --card=1 \
              -k \
              --name=$label \
              -o . \
              --freqlabel c1_f0 \
              --dir=$projectDir/../difx \
              --calconly \
              --startmjd=\$startmjd        
        fi
        """    
    
    stub:
        """
        mkdir c1_f0
        touch c1_f0/craftfrb0.im
        touch c1_f0/craftfrb0.calc
        """
}

process do_beamform {
    /*
        Produce a calibrated, beamformed fine spectrum for a particular
        antenna and polarisation from .vcraft voltages.

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
            imfile: path, calcfile: path
                The .im and .calc files that are used by craftcor_tab.py to 
                calculate geometric delays for beamforming
            pol: val
                One of "x" or "y" for the current polarisation being beamformed
            ant_idx: val
                Zero-based index of the antenna being beamformed
            flux_cal_solns: path
                Flux calibration solutions. These should be the same solutions 
                used to image the data and produce a position
            fcm: path
                fcm to use
            cand: val
                candidate file path, using val type so This can be reused with polcal, also file already
                exists before run, so don't need to wait for it
            dm: val
                DM of FRB

        Output
            pol, fine spectrum: tuple(val, path)
                The fine spectrum of the antenna polarisation

                The polarisation is included to be able to group outputs by 
                their polarisation
    */

    label 'celebi'

    input:
        val label
        val data
        tuple path(imfile), path(calcfile)
        each pol
        each ant_idx
        path flux_cal_solns
        path fcm
        val cand
        val dm

    output:
        tuple val(pol), path("${label}_frb_${ant_idx}_${pol}_f.npy"), emit: data
        env FFTLEN, emit: fftlen
		env bform_start_MJD, emit: bform_start_MJD
		
    script:
        """
        source /opt/setup_proc_container
        set -xu

        mkdir delays    # needed by craftcor_tab.py
        tar xvf $flux_cal_solns

        # High band FRBs need --uppersideband
        uppersideband=' '
        if [ "$params.uppersideband" = "true" ]; then
            uppersideband="--uppersideband"
        fi

        # Legacy compatibility: some very old FRBs need a hwfile
        hwfile=' '
        if [ ! "$params.hwfile" = "N/A" ]; then
            hwfile="--hwfile $params.hwfile"
        fi
        
        cropwins="-1.0"
        # Set cropping window        
        if [[ $label == "${params.label}" ]]; then
            if [[ "${params.longdata}" = "true" ]]; then
                cropwins="$params.frb_crop_width_s"
            fi
        else
            cropwins="$params.polcal_crop_width_s"
        fi
        
        candm=' '
        # Candidate file for cropping
        if [[ $label == "${params.label}" ]]; then
            candm="--snoopy=$cand --DM=$dm"
        fi

        echo "printing parameters"
        echo "\$hwfile \$candm"

        python3 $beamform_dir/craftcor_tab.py -d=$data \
            --parset=$fcm \
            --calcfile=$imfile \
            --aips_c bandpass*txt \
            --an=$ant_idx \
            --pol=$pol \
            -o=${label}_frb_${ant_idx}_${pol}_f.npy \
            -i=1 \
            --cpus=16 \
            --crop_width_s=\$cropwins \
            \$candm \$uppersideband \$hwfile

        rm TEMP*

        export FFTLEN=`cat fftlen`
        
        export bform_start_MJD=`cat corrected_start_MJD.txt`

        # if dirs do not exist, make them
        if [ ! -d ${params.out_dir}/htr ]; then
            mkdir ${params.out_dir}/htr
        fi

        if [ ! -d ${params.out_dir}/htr/info ]; then 
            mkdir ${params.out_dir}/htr/info
        fi

        # if txt file called ant_failed was produced, save it to htr/info folder of publish dir
        declare -a ffarr=(`find ./ -maxdepth 1 -name "*_failed.txt"`)
        if [[ \${#ffarr[@]} -gt 0 ]]; then
            cp *_failed.txt ${params.out_dir}/htr/info/.
        fi

        # only want to calculate this once, so choice of antenna id is arbitrary
        if [ $label == "${params.label}" ] && [ $ant_idx == 0 ]; then
            cp frb_crop_MJD.txt ${params.out_dir}/htr/info/frb_crop_MJD.txt
            cp corrected_start_MJD.txt ${params.out_dir}/htr/info/bform_start_MJD.txt
        fi

        if [ $ant_idx == 0 ]; then
            cp ant_vcraft_lengths.txt ${params.out_dir}/htr/info/${label}_ant_vcraft_lengths.txt
        fi

        # save txt file with cropping information of each antenna 
        cp ant_crop.txt ${params.out_dir}/htr/info/${label}_antcrop_${ant_idx}_${pol}_crop.txt
        """

    stub:
        """
        touch ${label}_frb_${ant_idx}_${pol}_f.npy
        export FFTLEN=100
        """
}

process sum_antennas {
    /*
        Sum fine spectra across antennas for a particular polarisation

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            pol, spectra: tuple(val, path)
                Polarisation and fine spectra files
        
        Output:
            pol, summed spectrum: tuple(val, path)
                Fine spectrum summed across all antennas in a single
                polarisation

                The polarisation is included to be able to group outputs by 
                their polarisation
    */

    label 'celebi'

    input:
        val label
        tuple val(pol), path(spectra)

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f.npy")

    script:
        """
        source /opt/setup_proc_container
        set -xu

        python3 $beamform_dir/sum.py \
                --f_dir . \
                -f=${label}_frb \
                -p=$pol \
                -o=${label}_frb_sum_${pol}_f.npy
        
        """
    
    stub:
        """
        touch ${label}_frb_sum_${pol}_f.npy
        """
}

process generate_deripple {
    /*
        Generate deripple coefficients based on number of samples in fine
        spectra

        Input
            fftlen: env
                Number of samples in fine spectra
        
        Output
            coeffs: path
                Derippling coefficients
    */

    label 'celebi'

    input:
        env FFTLEN
    
    output:
        path "deripple*npy", emit: coeffs

    script:
        """
        source /opt/setup_proc_container
        set -xu

        python3 $beamform_dir/generate_deripple.py \$FFTLEN $beamform_dir/.deripple_coeffs/ADE_R6_OSFIR.mat

        """
    
    stub:
        """
        touch deripple100.npy
        """
}

process deripple {
    /*
        Apply derippling coefficients to summed fine spectrum to cancel out
        systematic ripple        

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            pol, spectrum: tuple(val, path)
                Polarisation and fine spectrum file
            fftlen: env
                Number of samples in fine spectra
            coeffs: path
                Derippling coefficients
        
        Output:
            pol, derippled spectrum: tuple(val, path)
                Derippled fine spectrum summed across all antennas in a single
                polarisation

                The polarisation is included to be able to group outputs by 
                their polarisation
    */

    label 'celebi'

    input:
        val label
        tuple val(pol), path(spectrum)
        env FFTLEN
        path coeffs

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f_derippled.npy")

    script:
        """
        source /opt/setup_proc_container
        set -xu

        python3 $beamform_dir/deripple.py \
                -f=$spectrum \
                -l=\$FFTLEN \
                -o=${label}_frb_sum_${pol}_f_derippled.npy \
                --bw=$params.bw \
                -c=$coeffs \
                --cpus=1
        """
    
    stub:
        """
        touch ${label}_frb_sum_${pol}_f_derippled.npy
        """
}

process dedisperse {
    /*
        Coherently dedisperse a fine spectrum 

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            dm: val
                Dispersion measure to dedisperse to (pc/cm3)
            centre_freq: val
                Central frequency of fine spectrum (MHz)
            pol, spectrum: tuple(val, path)
                Polarisation and fine spectrum file
        
        Output:
            pol, dedispersed spectrum: tuple(val, path)
                Derippled, dedispersed fine spectrum summed across all antennas 
                in a single polarisation

                The polarisation is included to be able to group outputs by 
                their polarisation
    */

    label 'celebi'

    input:
        val label
        val dm
        val centre_freq
        tuple val(pol), path(spectrum)

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f_dedispersed_${dm}.npy")

    script:
        """
        source /opt/setup_proc_container 
        set -xu

        python3 $beamform_dir/dedisperse.py \
                -f=$spectrum \
                --DM=$dm \
                --f0=$centre_freq \
                --bw=$params.bw \
                -o=${label}_frb_sum_${pol}_f_dedispersed_${dm}.npy
        """
    
    stub:
        """
        touch ${label}_frb_sum_${pol}_f_dedispersed_${dm}.npy
        """
}

process ifft {
    /*
        Inverse fast Fourier transform fine spectrum to produce time series      

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            pol, spectrum: tuple(val, path)
                Polarisation and fine spectrum file
            dm: val
                Dispersion measure the data has been dedispersed to
        
        Output:
            pol_time_series: path
                ~3 ns dedispersed time series in a single polarisation    
    */

    label 'celebi'

    input:
        val label
        tuple val(pol), path(spectrum)
        val dm

    output:
        path("${label}_${pol}_t_${dm}.npy")

    script:
        """
        source /opt/setup_proc_container
        set -xu

        python3 $beamform_dir/ifft.py \
                -f=$spectrum \
                -o=${label}_${pol}_t_${dm}.npy

        # Copy the output into the publish_dir manually so Nextflow doesn't go over its
        # memory allocation
        if [ ! -d ${params.publish_dir}/${params.label}/htr ]; then
            mkdir ${params.publish_dir}/${params.label}/htr
        fi

        cp *_t_*.npy ${params.publish_dir}/${params.label}/htr/
        """

    stub:
        """
        touch ${label}_${pol}_t_${dm}.npy
        """
}

process generate_dynspecs {
    /*
        Generate Stokes parameter time series and dynamic spectra. 
        
        Generated time series will have (1/336 MHz) ~ 3 ns time resolution.
        Generated dynamic spectra will have 336 1 MHz channels at 1 us time
        resolution.

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            pol_time_series: path
                Two ~3 ns, dedispersed time series, one in each linear 
                polarisation
            centre_freq: val
                Central frequency of spectra (MHz)
            dm: val
                Dispersion measure the data has been dedispersed to
            pol_cal_solns: path
                Polarisation calibration solutions to be applied. If this is
                an empty file, polarisation calibration will not be applied.
        
        Output:
            data: path
                All .npy files created containing output Stokes parameter data
            dynspec_fnames: path
                File containing file names of dynamic spectra created
    */
    publishDir "${params.out_dir}/htr", mode: "copy"
    cpus 16
    label 'celebi'

    input:
        val label
        path pol_time_series
        val centre_freq
        val dm

    output:
        path "*dynspec*.npy", emit: data
        path "*fnames.txt", emit: dynspec_fnames
        path "*.png"
        path "*.npy"

    script:
        """
        source /opt/setup_proc_container
        set -xu
		
		touch dummy.png
		
		bfautoflg=' '
		if [ "$params.bform_autoflag" == "true" ]; then
            bfautoflg="--do_chanflag"
        fi
        
        do_fbline_corr=' '   
        if [ "$params.frb_bline" == "true" ]; then
            do_fbline_corr="--bline"
        fi  
        
        do_pbline_corr=' '   
        if [ "$params.pcal_bline" == "true" ]; then
            do_pbline_corr="--bline"
        fi   
		
        if [[ $label == "${params.label}_polcal" ]]; then
            type="polcal"
            # this feels illegal
            MJD1=\$(echo \$(<$params.snoopy) | cut -d ' ' -f 21)

            python3 $beamform_dir/make_dynspec.py \
                    -x=${label}_X_t_${dm}.npy \
                    -y=${label}_Y_t_${dm}.npy \
                    \$do_pbline_corr \
                    --ofile=${label}_@_dynspec_${dm}.npy \
                    --chanlists=$projectDir/../flagging \
                    \$bfautoflg \
                    --pulsar \
                    --MJD0=$params.polcal_MJD0 \
                    --MJD1=\$MJD1 \
                    --F0=$params.polcal_F0 \
                    --F1=$params.polcal_F1 \
                    --DM=$dm \
                    --cfreq=$centre_freq \
                    --bw=$params.bw \
                    --sigma=$params.polcal_dynspec_sigma \
                    --baseline=$params.polcal_baseline \
                    --tN=$params.polcal_dynspec_tN \
                    --guard=$params.polcal_dynspec_guard  
        else
            type="frb"
            python3 $beamform_dir/make_dynspec.py \
                    -x=${label}_X_t_${dm}.npy \
                    -y=${label}_Y_t_${dm}.npy \
                    \$do_pbline_corr \
                    --ofile=${label}_@_dynspec_${dm}.npy \
                    --chanlists=$projectDir/../flagging \
                    \$bfautoflg \
                    --sigma=$params.frb_dynspec_sigma \
                    --baseline=$params.frb_baseline \
                    --tN=$params.frb_dynspec_tN \
                    --guard=$params.frb_dynspec_guard 
        fi
		
		# copy flagging files
        cp flagged_channels.npy ${label}_flagged_channels.npy
		
        echo "${label}_I_dynspec_${dm}.npy" > dynspec_fnames.txt
        echo "${label}_Q_dynspec_${dm}.npy" >> dynspec_fnames.txt
        echo "${label}_U_dynspec_${dm}.npy" >> dynspec_fnames.txt
        echo "${label}_V_dynspec_${dm}.npy" >> dynspec_fnames.txt
        
        # saving files
        if [ -f fail_bline.png ]; then
            cp fail_bline.png ${params.out_dir}/htr/${label}_fail_bline.png
            cp fail_bline.npy ${params.out_dir}/htr/${label}_fail_bline.npy
            echo "BASELINE CORRECTION FAILED for ${label}, CHECK ${label}_fail_bline files in htr directory"
        fi
        """
    
    stub:
        """
        touch stub.npy
        touch stub.txt
        """
}


workflow beamform {
    /*
        Workflow to produce beamformed voltage data from vcraft voltages.

        Take
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
            pos: path
                File containing JMFIT statistics of position to beamform on
            flux_cal_solns: path
                Flux calibration solutions. These should be the same solutions 
                used to image the data and produce a position
            pol_cal_solns: path
                Polarisation calibration solutions to be applied. If this is
                an empty file, polarisation calibration will not be applied.
            dm: val
                Dispersion measure to dedisperse to (pc/cm3)
            centre_freq: val
                Central frequency of fine spectrum (MHz)
            ds_args: val
                String containing arguments to be passed to dynspecs.py. Use
                this to specify which Stokes parameters and data types (time
                series or dynamic spectrum) to generate.
            nants: val
                Number of antennas available in the data
            fcm: path
                fcm to use
            cand: val
                candidate file path, using val type so This can be reused with polcal, also file already
                exists before run, so don't need to wait for it
        
        Emit
        	xy: 
        		beamformed X,Y voltages
        	pre_dedisp:
        		Outputs of derippling
        	bform_start_MJD:
        		Start timestamp for beamformed data
              
    */
    take:
        label               // FRB label
        data                // data 
        pos                 // frb position
        flux_cal_solns      // flux cal solutions
        pol_cal_solns       // pol cal solutions
        dm                  // DM
        centre_freq         // central frequency
        nants               // number of antennas
        fcm                 // fcm file
        cand                // path to cand file
    
    main:
        // preliminaries
        calcfiles = create_calcfiles(label, data, pos, fcm)

        antennas = Channel
            .of(0..nants-1)

        //antennas.view()

        // processing
        
        // apply delays and calibration solutions to each antenna/pol fine spectra, align each antenna
        do_beamform(
            label, data, calcfiles, polarisations, antennas, flux_cal_solns, fcm, cand, dm
        )

        // filter antenna to make sure non-empty data is being beamformed
        filter_antenna(label, do_beamform.out.data)

        // sum filtered antenna data together to get summed X and Y polarisation fine spectra -> beamforming
        sum_antennas(label, filter_antenna.out.filtered_ant.groupTuple())

        // calculate derriple coefficients
        coeffs = generate_deripple(do_beamform.out.fftlen.first())

        // apply deripple coefficients
        deripple(label, sum_antennas.out, do_beamform.out.fftlen.first(), coeffs)

        // coherently dedisperse fine spectra
        dedisperse(label, dm, centre_freq, deripple.out)

        // inverse FFT back to complex time series data
        ifft(label, dedisperse.out, dm)
        xy = ifft.out.collect()

        // if FRB, apply polcal solutions to x and y data products
        if ((label == "${params.label}") && !params.nopolcal) {
            xy = apply_pol_cal_solns(label, xy, pol_cal_solns, centre_freq, dm).calib_data
            label="${params.label}_calib"
        }
    
    emit:
        xy
        pre_dedisp = deripple.out
        bform_start_MJD = do_beamform.out.bform_start_MJD.first()
}


workflow gen_dspec {
    /*
        Workflow to produce high-time resolution time series and dynamic
        spectra across Stokes IQUV 

        Take
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            xy: val
                Path to XY beamformed votages
            dm: val
                Dispersion measure to dedisperse to (pc/cm3)
            centre_freq: val
                Central frequency of fine spectrum (MHz)
        
        Emit
        	dynamic_fnames:
        		Filenames for dynamic spectra
            htr_data: path
                Numpy files containing Stokes time series and dynamic spectra    
    */
    take:
        label               // FRB label
        xy                  // XY voltages        
        dm                  // DM
        centre_freq         // central frequency
    
    main:        
        // generate stokes I, Q, U and V dynamic spectra
        if ((label == "${params.label}") && !params.nopolcal) {
            alabel="${params.label}_calib"
        }
        else {
            alabel=label
        }
        generate_dynspecs(alabel, xy, centre_freq, dm)
    
    emit:
        dynspec_fnames = generate_dynspecs.out.dynspec_fnames
        htr_data = generate_dynspecs.out.data
}
