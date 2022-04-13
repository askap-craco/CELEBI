nextflow.enable.dsl=2

include { get_num_ants } from './utils'
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
        export CRAFTCATDIR="."

        ra=\$(grep "Actual RA" $pos)
        ra=\${ra:22}
        dec=\$(grep "Actual Dec" $pos)
        dec=\${dec:22}

        # Run processTimeStep.py with the --calconly flag to stop once calcfile is written
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

        echo "$localise_dir/processTimeStep.py \$args"
        $localise_dir/processTimeStep.py \$args
        """    
}

process do_beamform {
    cpus 4  // need to limit max number of parallel processes to 4 for memory

    input:
        val label
        val data
        tuple path(imfile), path(calcfile)
        tuple val(pol), val(ant_idx)
        path flux_cal_solns
        val num_ints
        val int_len
        val offset
        val fcm

    output:
        tuple val(pol), path("${label}_frb_${ant_idx}_${pol}_f.npy")

    script:
        """
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

        # Legacy compatibility: some very old FRBs need a hwfile
        if [ ! "$params.hwfile" = "N/A" ]; then
            args="\$args --hwfile $params.hwfile"
        fi

        python $beamform_dir/craftcor_tab.py \$args
        rm TEMP*
        """
}

process sum {
    input:
        val label
        tuple val(pol), path(spectra)

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f.npy")

    script:
        """
        args="--f_dir ."
        args="\$args -f ${label}_frb"
        args="\$args -p $pol"
        args="\$args -o ${label}_frb_sum_${pol}_f.npy"

        python3 $beamform_dir/sum.py \$args
        """
}

process deripple {
    cpus 1
    memory '60 GB'

    input:
        val label
        val int_len
        tuple val(pol), path(spectrum)

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f_derippled.npy")

    """
    module load python/2.7.14
    module load numpy/1.16.3-python-2.7.14
    module load scipy/1.0.0-python-2.7.14

    fftlen=\$(( $int_len * 64 ))

    if [ ! -d $beamform_dir/.deripple_coeffs ]; then
        mkdir $beamform_dir/.deripple_coeffs
    fi

    args="-f $spectrum"
    args="\$args -l \$fftlen"
    args="\$args -o ${label}_frb_sum_${pol}_f_derippled.npy"
    args="\$args -c $beamform_dir/.deripple_coeffs"

    python $beamform_dir/deripple.py \$args
    """
}

process dedisperse {
    executor 'slurm'
    cpus 1
    time '5m'
    memory '64 GB'

    input:
        val label
        val dm
        val centre_freq
        tuple val(pol), path(spectrum)

    output:
        tuple val(pol), path("${label}_frb_sum_${pol}_f_dedispersed.npy")

    script:
        """
        module load gcc/7.3.0
        module load openmpi/3.0.0
        module load python/3.7.4
        module load numpy/1.18.2-python-3.7.4

        args="-f $spectrum"
        args="\$args --DM $DM"
        args="\$args --f0 $centre_freq"
        args="\$args --bw 336"
        args="\$args -o ${label}_frb_sum_${pol}_f_dedispersed.npy"

        python3 $beamform_dir/dedisperse.py \$args
        """
}

process ifft {
    input:
    val label
    tuple val(pol), path(spectrum)
    path pol_cal_solns

    output:
    path("${label}_frb_sum_${pol}_t.npy")

    """
    args="-f $spectrum"
    args="\$args -o ${label}_frb_sum_${pol}_t.npy"

    python3 $beamform_dir/ifft.py \$args
    """
}

process generate_dynspecs {
    executor 'slurm'
    cpus 1
    time '1h'
    memory '64 GB'

    input:
        val label
        path pol_time_series

    output:
        path "${label}_frb_fulltimeres.tar.gz"

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    args="-x ${label}_frb_sum_x_t.npy"
    args="\$args -y ${label}_frb_sum_y_t.npy"
    args="\$args -o ${label}_frb_sum_!_@.npy"

    python3 $beamform_dir/dynspecs.py \$args

    tar -czvhf ${label}_frb_fulltimeres.tar.gz *.npy
    """
}

workflow beamform {
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
    
    main:
        // preliminaries
        startmjd = get_startmjd(data)
        calcfiles = create_calcfiles(label, data, startmjd, pos, fcm)

        // processing
        do_beamform(
            label, data, calcfiles, polarisations.combine(antennas), flux_cal_solns,
            num_ints, int_len, offset, fcm
        )
        sum(label, do_beamform.out.groupTuple())
        deripple(label, int_len, sum.out)
        dedisperse(label, dm, centre_freq, deripple.out)
        ifft(label, dedisperse.out, pol_cal_solns)
        generate_dynspecs(label, ifft.out.collect())
    
    emit:
        htr_data = generate_dynspecs.out
}
