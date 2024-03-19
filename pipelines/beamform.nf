/*
    Processes and workflow related to beamforming and producing high-time 
    resolution data products
*/

nextflow.enable.dsl=2

include { get_start_mjd } from './correlate'
include { apply_pol_cal_solns } from './calibration'

params.pols = ['X', 'Y']
polarisations = Channel
    .fromList(params.pols)

params.hwfile = "N/A"

localise_dir = "$baseDir/../localise/"
beamform_dir = "$baseDir/../beamform/"
params.uppersideband = false
params.out_dir = "${params.publish_dir}/${params.label}"

params.bw = 336 /*Default value*/


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
    input:
        val label
        val data
        path pos
        path fcm

    output:
        tuple path("c1_f0/craftfrb*.im"), path("c1_f0/craftfrb*.calc")

    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        #fi
        ml apptainer
        set -a
        set -o allexport
        startmjd=`python3 $localise_dir/get_start_mjd.py $data`

        export CRAFTCATDIR="."

        ra=\$(grep "Actual RA" $pos)
        ra=\${ra:22}
        dec=\$(grep "Actual Dec" $pos)
        dec=\${dec:22}

        # Run processTimeStep.py with the --calconly flag to stop once calcfile is 
        # written
        args="-t $data"
        args="\$args --ra \$ra"
        args="\$args -d\$dec"
        args="\$args -f $fcm"
        args="\$args -b 4"
        args="\$args --card 1"
        args="\$args -k"
        args="\$args --name=$label"
        args="\$args -o ."
        args="\$args --freqlabel c1_f0"
        args="\$args --dir=$baseDir/../difx"
        args="\$args --calconly"
        args="\$args --startmjd \$startmjd"

        echo "python3 $localise_dir/processTimeStep.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/processTimeStep.py \$args'
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

        Output
            pol, fine spectrum: tuple(val, path)
                The fine spectrum of the antenna polarisation

                The polarisation is included to be able to group outputs by 
                their polarisation
    */
    input:
        val label
        val data
        tuple path(imfile), path(calcfile)
        each pol
        each ant_idx
        path flux_cal_solns
        path fcm

    output:
        tuple val(pol), path("${label}_frb_${ant_idx}_${pol}_f.npy"), emit: data
        env FFTLEN, emit: fftlen

    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
            #. $launchDir/../setup_beamform
        #    . $launchDir/../setup_proc
        #fi
        mkdir delays    # needed by craftcor_tab.py
        tar xvf $flux_cal_solns
        ml apptainer
        set -a
        set -o allexport
        args="-d $data"
        args="\$args --parset $fcm"
        args="\$args --calcfile $imfile"
        args="\$args --aips_c bandpass*txt"
        args="\$args --an $ant_idx"
        args="\$args --pol $pol"
        args="\$args -o ${label}_frb_${ant_idx}_${pol}_f.npy"
        args="\$args -i 1"
        args="\$args --cpus=16"

        # High band FRBs need --uppersideband
        if [ "$params.uppersideband" = "true" ]; then
            args="\$args --uppersideband"
        fi

        # Legacy compatibility: some very old FRBs need a hwfile
        if [ ! "$params.hwfile" = "N/A" ]; then
            args="\$args --hwfile $params.hwfile"
        fi

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/craftcor_tab.py \$args'
        rm TEMP*

        export FFTLEN=`cat fftlen`
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
    input:
        val label
        tuple val(pol), path(spectra)

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f.npy")

    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
            #. $launchDir/../setup_beamform
        #    . $launchDir/../setup_proc
        #fi
        ml apptainer
        set -a
        set -o allexport
        args="--f_dir ."
        args="\$args -f ${label}_frb"
        args="\$args -p $pol"
        args="\$args -o ${label}_frb_sum_${pol}_f.npy"

        echo "python3 $beamform_dir/sum.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/sum.py \$args'
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
    input:
        env FFTLEN
    
    output:
        path "deripple*npy", emit: coeffs

    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
        #    module load gcc/9.2.0
        #    module load openmpi/4.0.2
        #    module load python/3.7.4
        #    module load numpy/1.18.2-python-3.7.4
        #    module load scipy/1.6.0-python-3.7.4
        #    export PYTHONPATH=\$PYTHONPATH:/fred/oz002/askap/craft/craco/python/lib/python3.7/site-packages/
        #fi
        ml apptainer
        set -a
        set -o allexport
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/generate_deripple.py \$FFTLEN $beamform_dir/.deripple_coeffs/ADE_R6_OSFIR.mat'

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
    input:
        val label
        tuple val(pol), path(spectrum)
        env FFTLEN
        path coeffs

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f_derippled.npy")

    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
        #    module load gcc/9.2.0
        #    module load openmpi/4.0.2
        #    module load python/3.7.4
        #    module load numpy/1.18.2-python-3.7.4
        #    module load scipy/1.6.0-python-3.7.4
        #    module load joblib/0.11
        #    export PYTHONPATH=\$PYTHONPATH:/fred/oz002/askap/craft/craco/python/lib/python3.7/site-packages/
        #fi
        ml apptainer
        set -a
        set -o allexport
        args="-f $spectrum"
        args="\$args -l \$FFTLEN"
        args="\$args -o ${label}_frb_sum_${pol}_f_derippled.npy"
        args="\$args --bw $params.bw"
        args="\$args -c $coeffs"
        args="\$args --cpus 1"

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/deripple.py \$args'
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

    input:
        val label
        val dm
        val centre_freq
        tuple val(pol), path(spectrum)

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f_dedispersed_${dm}.npy")

    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
            #. $launchDir/../setup_beamform
        #    . $launchDir/../setup_proc
        #fi
        ml apptainer
        set -a
        set -o allexport
        args="-f $spectrum"
        args="\$args --DM $dm"
        args="\$args --f0 $centre_freq"
        args="\$args --bw $params.bw"
        args="\$args -o ${label}_frb_sum_${pol}_f_dedispersed_${dm}.npy"

        echo "python3 $beamform_dir/dedisperse.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/dedisperse.py \$args'
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
    input:
        val label
        tuple val(pol), path(spectrum)
        val dm

    output:
        path("${label}_${pol}_t_${dm}.npy")

    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
        #    module load ipp/2018.2.199
        #    module load gcc/9.2.0
        #    module load openmpi/4.0.2
        #    module load gsl/2.5
        #    module load python/3.7.4
        #    module load numpy/1.18.2-python-3.7.4
        #    module load scipy/1.4.1-python-3.7.4
        #fi
        ml apptainer
        set -a
        set -o allexport

        args="-f $spectrum"
        args="\$args -o ${label}_${pol}_t_${dm}.npy"

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/ifft.py \$args'

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

    input:
        val label
        path pol_time_series
        val centre_freq
        val dm

    output:
        path "*.npy", emit: data
        path "*.txt", emit: dynspec_fnames

    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
            #. $launchDir/../setup_beamform
        #    . $launchDir/../setup_proc
        #fi
        ml apptainer
        set -a
        set -o allexport

        args="-x ${label}_X_t_${dm}.npy"
        args="\$args -y ${label}_Y_t_${dm}.npy"
        args="\$args --bline"
        args="\$args --ofile ${label}_@_dynspec_${dm}.npy"

        if [[ $label == "${params.label}_polcal" ]]; then
            type="polcal"
        else
            type="frb"
        fi
        

        if [ \$type == "polcal" ]; then 
            args="\$args --pulsar"
            args="\$args --MJD0 $params.polcal_MJD0"

            # this feels illegal
            MJD1=\$(echo \$(<$params.snoopy) | cut -d ' ' -f 21)
            echo \$MJD1
            args="\$args --MJD1 \$MJD1"

            args="\$args --F0 $params.polcal_F0"
            args="\$args --F1 $params.polcal_F1"
            args="\$args --DM $dm"
            args="\$args --cfreq $centre_freq"
            args="\$args --bw $params.bw"
            args="\$args --sigma $params.polcal_dynspec_sigma"
            args="\$args --baseline $params.polcal_baseline"
            args="\$args --tN $params.polcal_dynspec_tN"
            args="\$args --guard $params.polcal_dynspec_guard"
        
        elif [ \$type == "frb" ]; then
            args="\$args --sigma $params.frb_dynspec_sigma"
            args="\$args --baseline $params.frb_baseline"
            args="\$args --tN $params.frb_dynspec_tN"
            args="\$args --guard $params.frb_dynspec_guard"
        fi

        echo "python3 $beamform_dir/make_dynspec.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/make_dynspec.py \$args'


        echo "${label}_I_dynspec_${dm}.npy" > dynspec_fnames.txt
        echo "${label}_Q_dynspec_${dm}.npy" >> dynspec_fnames.txt
        echo "${label}_U_dynspec_${dm}.npy" >> dynspec_fnames.txt
        echo "${label}_V_dynspec_${dm}.npy" >> dynspec_fnames.txt

        """
    
    stub:
        """
        touch stub.npy
        touch stub.txt
        """
}


workflow beamform {
    /*
        Workflow to produce high-time resolution time series and dynamic
        spectra across Stokes IQUV from vcraft voltages.

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
        
        Emit
            htr_data: path
                Numpy files containing Stokes time series and dynamic spectra    
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
    
    main:
        // preliminaries
        calcfiles = create_calcfiles(label, data, pos, fcm)

        antennas = Channel
            .of(0..nants-1)

        // processing
        do_beamform(
            label, data, calcfiles, polarisations, antennas, flux_cal_solns, fcm
        )
        sum_antennas(label, do_beamform.out.data.groupTuple())
        coeffs = generate_deripple(do_beamform.out.fftlen.first())
        deripple(label, sum_antennas.out, do_beamform.out.fftlen.first(), coeffs)
        dedisperse(label, dm, centre_freq, deripple.out)
        ifft(label, dedisperse.out, dm)
        xy = ifft.out.collect()

        // if FRB, apply polcal solutions
        if (label == "${params.label}") {
            xy = apply_pol_cal_solns(label, xy, pol_cal_solns, centre_freq, dm).calib_data
            label="${params.label}_calib"
        }

        generate_dynspecs(label, xy, centre_freq, dm)
    
    emit:
        dynspec_fnames = generate_dynspecs.out.dynspec_fnames
        htr_data = generate_dynspecs.out.data
        xy
        pre_dedisp = deripple.out
}
