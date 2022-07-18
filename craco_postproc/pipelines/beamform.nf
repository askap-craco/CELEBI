/*
    Processes and workflow related to beamforming and producing high-time 
    resolution data products
*/

nextflow.enable.dsl=2

include { get_startmjd } from './correlate'

params.pols = ['x', 'y']
polarisations = Channel
    .fromList(params.pols)
params.nants = 2
antennas = Channel
    .of(0..params.nants-1)

params.hwfile = "N/A"

localise_dir = "$baseDir/../localise/"
beamform_dir = "$baseDir/../beamform/"

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
            startmjd: val
                Earliest data start time in Modified Julian Day (MJD)
            pos: path
                File containing JMFIT statistics of position to beamform on
            fcm: val
                Absolute path to fcm (hardware delays) file
        
        Output
            Delay files: tuple(path, path)
                The .im and .calc files that are used by craftcor_tab.py to 
                calculate geometric delays for beamforming
    */
    input:
        val label
        val data
        val startmjd
        path pos
        val fcm

    output:
        tuple path("c1_f0/craftfrb*.im"), path("c1_f0/craftfrb*.calc")

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi

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
        args="\$args --startmjd $startmjd"

        echo "python3 $localise_dir/processTimeStep.py \$args"
        python3 $localise_dir/processTimeStep.py \$args
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
            num_ints: val
                Number of integrations
            int_len: val
                Integration length (units unclear?)
            offset: val
                Integration start offset (units unclear, not same as int_len)
            fcm: val
                Absolute path to fcm (hardware delays) file

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
        val num_ints
        val int_len
        val offset
        val fcm

    output:
        tuple val(pol), path("${label}_frb_${ant_idx}_${pol}_f.npy"), emit: data
        env FFTLEN, emit: fftlen

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_beamform
        fi
        mkdir delays    # needed by craftcor_tab.py
        tar xvf $flux_cal_solns

        args="-i $num_ints"
        args="\$args -n $int_len"
        args="\$args --offset $offset"
        args="\$args -d $data"
        args="\$args --parset $fcm"
        args="\$args --calcfile $imfile"
        args="\$args --aips_c bandpass*txt"
        args="\$args --an $ant_idx"
        args="\$args --pol $pol"
        args="\$args -o ${label}_frb_${ant_idx}_${pol}_f.npy"
        args="\$args --tab"
        args="\$args --cpus=16"

        # Legacy compatibility: some very old FRBs need a hwfile
        if [ ! "$params.hwfile" = "N/A" ]; then
            args="\$args --hwfile $params.hwfile"
        fi

        python3 $beamform_dir/craftcor_tab.py \$args
        rm TEMP*

        export FFTLEN=`cat fftlen`
        """
}

process sum {
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
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_beamform
        fi
        args="--f_dir ."
        args="\$args -f ${label}_frb"
        args="\$args -p $pol"
        args="\$args -o ${label}_frb_sum_${pol}_f.npy"

        python3 $beamform_dir/sum.py \$args
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
        path "*npy", emit: coeffs

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            module load gcc/9.2.0
            module load openmpi/4.0.2
            module load python/3.7.4
            module load numpy/1.18.2-python-3.7.4
            module load scipy/1.6.0-python-3.7.4
            export PYTHONPATH=\$PYTHONPATH:/fred/oz002/askap/craft/craco/python/lib/python3.7/site-packages/
        fi

        python3 $beamform_dir/generate_deripple.py \$FFTLEN $beamform_dir/.deripple_coeffs/ADE_R6_OSFIR.mat
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
            int_len: val
                Integration length (units unclear?)
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
        val int_len
        tuple val(pol), path(spectrum)
        env FFTLEN
        path coeffs

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f_derippled.npy")

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            module load gcc/9.2.0
            module load openmpi/4.0.2
            module load python/3.7.4
            module load numpy/1.18.2-python-3.7.4
            export PYTHONPATH=\$PYTHONPATH:/fred/oz002/askap/craft/craco/python/lib/python3.7/site-packages/
        fi

        args="-f $spectrum"
        args="\$args -l \$FFTLEN"
        args="\$args -o ${label}_frb_sum_${pol}_f_derippled.npy"
        args="\$args -c $coeffs"
        args="\$args --cpus 1"

        python3 $beamform_dir/deripple.py \$args
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
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_beamform
        fi
        args="-f $spectrum"
        args="\$args --DM $dm"
        args="\$args --f0 $centre_freq"
        args="\$args --bw 336"
        args="\$args -o ${label}_frb_sum_${pol}_f_dedispersed_${dm}.npy"

        echo "python3 $beamform_dir/dedisperse.py \$args"
        python3 $beamform_dir/dedisperse.py \$args
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
            pol_cal_solns: path
                Polarisation calibration solutions to be applied. If this is
                an empty file, polarisation calibration will not be applied.
            dm: val
                Dispersion measure the data has been dedispersed to
        
        Output:
            pol_time_series: path
                ~3 ns dedispersed time series in a single polarisation    
    */
    input:
        val label
        tuple val(pol), path(spectrum)
        path pol_cal_solns
        val dm

    output:
        path("${label}_frb_sum_${pol}_t_${dm}.npy")

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            module load ipp/2018.2.199
            module load gcc/9.2.0
            module load openmpi/4.0.2
            module load gsl/2.5
            module load python/3.7.4
            module load numpy/1.18.2-python-3.7.4
            module load scipy/1.4.1-python-3.7.4
        fi
        args="-f $spectrum"
        args="\$args -o ${label}_frb_sum_${pol}_t_${dm}.npy"

        python3 $beamform_dir/ifft.py \$args

        # Copy the output into the publish_dir manually so Nextflow doesn't go over its
        # memory allocation
        if [ ! -d ${params.publish_dir}/${params.label}/htr ]; then
            mkdir ${params.publish_dir}/${params.label}/htr
        fi

        cp *_t_*.npy ${params.publish_dir}/${params.label}/htr/
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
            ds_args: val
                String containing arguments to be passed to dynspecs.py. Use
                this to specify which Stokes parameters and data types (time
                series or dynamic spectrum) to generate.
            dm: val
                Dispersion measure the data has been dedispersed to
        
        Output:
            data: path
                All .npy files created containing output Stokes parameter data
            dynspec_fnames: path
                File containing file names of dynamic spectra created
    */
    publishDir "${params.publish_dir}/${params.label}/htr", mode: "copy"
    cpus 16

    input:
        val label
        path pol_time_series
        val ds_args
        val dm

    output:
        path "*.npy", emit: data
        path "*.txt", emit: dynspec_fnames

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_beamform
        fi
        args="-x ${label}_frb_sum_x_t_${dm}.npy"
        args="\$args -y ${label}_frb_sum_y_t_${dm}.npy"
        args="\$args -o ${label}_frb_sum_!_@_${dm}.npy" 

        echo "python3 $beamform_dir/dynspecs.py \$args $ds_args"
        python3 $beamform_dir/dynspecs.py \$args $ds_args
        """
}

process plot {
    /*
        Plot dynamic spectra across different time resolutions to produce
        summary plot.

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            fnames_file: path
                File containing file names of dynamic spectra
            dynspecs: path
                Stokes parameter dynamic spectra
            centre_freq: val
                Central frequency of fine spectrum (MHz)
            dm: val
                Dispersion measure the data has been dedispersed to
        
        Output:
            plot: path
                Plotted dynamic spectra across different time resolutions
    */
    publishDir "${params.publish_dir}/${params.label}/htr", mode: "copy"

    input:
        val label
        path fnames_file
        path dynspecs
        val centre_freq
        val dm
    
    output:
        path "*.png"
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_beamform
        fi
        args="-s $fnames_file"
        args="\$args -f $centre_freq"
        args="\$args -l $label"
        args="\$args -d $dm"

        python3 $beamform_dir/plot.py \$args
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
            fcm: val
                Absolute path to fcm (hardware delays) file
            pos: path
                File containing JMFIT statistics of position to beamform on
            flux_cal_solns: path
                Flux calibration solutions. These should be the same solutions 
                used to image the data and produce a position
            pol_cal_solns: path
                Polarisation calibration solutions to be applied. If this is
                an empty file, polarisation calibration will not be applied.
            num_ints: val
                Number of integrations
            int_len: val
                Integration length (units unclear?)
            offset: val
                Integration start offset (units unclear, not same as int_len)
            dm: val
                Dispersion measure to dedisperse to (pc/cm3)
            centre_freq: val
                Central frequency of fine spectrum (MHz)
            ds_args: val
                String containing arguments to be passed to dynspecs.py. Use
                this to specify which Stokes parameters and data types (time
                series or dynamic spectrum) to generate.
        
        Emit
            htr_data: path
                Numpy files containing Stokes time series and dynamic spectra    
    */
    take:
        label   // val
        data    // val
        fcm // val
        pos // path
        flux_cal_solns  // path
        pol_cal_solns   // path
        num_ints    // val
        int_len // val
        offset  // val
        dm  // val
        centre_freq // val
        ds_args // val
    
    main:
        // preliminaries
        startmjd = get_startmjd(data)
        calcfiles = create_calcfiles(label, data, startmjd, pos, fcm)

        // processing
        do_beamform(
            label, data, calcfiles, polarisations, antennas, flux_cal_solns,
            num_ints, int_len, offset, fcm
        )
        sum(label, do_beamform.out.data.groupTuple())
        coeffs = generate_deripple(do_beamform.out.fftlen.first())
        deripple(label, int_len, sum.out, do_beamform.out.fftlen.first(), coeffs)
        dedisperse(label, dm, centre_freq, deripple.out)
        ifft(label, dedisperse.out, pol_cal_solns, dm)
        generate_dynspecs(label, ifft.out.collect(), ds_args, dm)
        plot(label, generate_dynspecs.out.dynspec_fnames, generate_dynspecs.out.data, centre_freq, dm)
    
    emit:
        htr_data = generate_dynspecs.out.data
}
