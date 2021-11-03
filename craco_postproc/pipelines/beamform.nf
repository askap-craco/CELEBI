#!/usr/bin/env nextflow

params.user = "dscott"

params.cal_file = "$baseDir/calibrators.dat"

params.skip_polcal = false

antennas_polcal = Channel
    .from(0..params.nant_polcal-1)

antennas_frb = Channel
    .from(0..params.nant_frb-1)

params.pols = ['x', 'y']
polarisations = Channel
    .fromList(params.pols)
    .into{polarisations_polcal; polarisations_frb}

sourceDir="$baseDir/beamform"

process create_calc_polcal {
    output:
    tuple path("c1_f0/craftfrb*.im"), path("c1_f0/craftfrb*.calc") into calcfiles_polcal

    when:
    params.skip_polcal == false

    """
    source $baseDir/setup_difx $params.user

    export CRAFTCATDIR="."

    # Run processTimeStep.py with the --calconly flag to stop once calcfile is written
    args="-t $params.data_polcal"
    args="\$args --ra `grep $params.polcal_name $params.cal_file | cut -d " " -f 2`"
    args="\$args -d`grep $params.polcal_name $params.cal_file | cut -d " " -f 3`"
    args="\$args -f $params.fcm"
    args="\$args -b 4"
    args="\$args --card 1"
    args="\$args -k"
    args="\$args --slurm"
    args="\$args --name=${params.label}_polcal"
    args="\$args -o ."
    args="\$args --freqlabel c1_f0"
    args="\$args --dir=$baseDir/difx"
    args="\$args --calconly"

    $baseDir/craftpy2/processTimeStep.py \$args
    """
}

process get_startmjd_frb {
    output:
    env startmjd into startmjd_frb

    """
    module load python/3.8.5

    startmjd=`python3 $baseDir/scripts/get_start_mjd.py $params.data_frb`
    echo "\$startmjd"
    """
}

process create_calc_frb {
    input:
    env startmjd from startmjd_frb

    output:
    tuple path("c1_f0/craftfrb*.im"), path("c1_f0/craftfrb*.calc") into created_calcfiles_frb

    when:
    params.calcdir == ""

    """
    source $baseDir/setup_difx $params.user

    export CRAFTCATDIR="."

    # Run processTimeStep.py with the --calconly flag to stop once calcfile is written
    args="-t $params.data_frb"
    args="\$args --ra $params.ra_frb"
    args="\$args -d$params.dec_frb"
    args="\$args -f $params.fcm"
    args="\$args -b 4"
    args="\$args --card 1"
    args="\$args -k"
    args="\$args --slurm"
    args="\$args --name=${params.label}_frb"
    args="\$args -o ."
    args="\$args --freqlabel c1_f0"
    args="\$args --dir=$baseDir/difx"
    args="\$args --calconly"
    args="\$args --startmjd \$startmjd"

    $baseDir/craftpy2/processTimeStep.py \$args
    """
}

process get_calc_frb {
    output:
    tuple path("craftfrb*.im"), path("craftfrb*.calc") into got_calcfiles_frb

    when:
    params.calcdir != ""

    """
    cp $params.calcdir/*.calc .
    cp $params.calcdir/*.im .
    """
}

process beamform_polcal {
    executor 'slurm'
    cpus 4
    time '20m'
    memory '32 GB'

    input:
    tuple path(imfile), path(calcfile) from calcfiles_polcal
    each pol from polarisations_polcal
    each antnum from antennas_polcal

    output:
    tuple val(pol), path("${params.label}_polcal_${antnum}_${pol}_f.npy") into spectra_polcal

    """
    module load python/2.7.14 
    module load numpy/1.16.3-python-2.7.14 
    module load scipy/1.0.0-python-2.7.14 
    module load astropy/2.0.3-python-2.7.14 
    module load matplotlib/2.2.2-python-2.7.14 
    module load joblib/0.11

    if [ ! -d delays ]; then
        mkdir delays    # needed by craftcor_tab.py
    fi

    args="-i $params.numints_polcal"
    args="\$args -n $params.intlen_polcal"
    args="\$args --offset $params.offset_polcal"
    args="\$args -d $params.data_polcal"
    args="\$args --parset $params.fcm"
    args="\$args --calcfile $imfile"
    args="\$args --aips_c $params.bandpass"
    args="\$args --an $antnum"
    args="\$args --pol $pol"
    args="\$args -o ${params.label}_polcal_${antnum}_${pol}_f.npy"
    args="\$args --tab"
    if [ ! "${params.hwfile}" = "" ]; then
        args="\$args --hwfile $params.hwfile"
    fi

    python $sourceDir/craftcor_tab.py \$args
    """
}

process sum_polcal {
    executor 'slurm'
    cpus 1
    time '15m'
    memory '16 GB'

    input:
    tuple val(pol), path(spectra) from spectra_polcal.groupTuple()

    output:
    tuple val(pol), path("${params.label}_polcal_sum_${pol}_f.npy") into summed_spectrum_polcal

    """
    module load gcc/7.3.0 
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    args="--f_dir ."
    args="\$args -f ${params.label}_polcal"
    args="\$args -p $pol"
    args="\$args -o ${params.label}_polcal_sum_${pol}_f.npy"

    python3 $sourceDir/sum.py \$args
    """
}

process deripple_polcal {
    executor 'slurm'
    cpus 1
    time '45m'
    memory '64 GB'

    input:
    tuple val(pol), path(spectrum) from summed_spectrum_polcal

    output:
    tuple val(pol), path("${params.label}_polcal_sum_${pol}_f_derippled.npy") into derippled_spectrum_polcal

    """
    module load python/2.7.14
    module load numpy/1.16.3-python-2.7.14
    module load scipy/1.0.0-python-2.7.14

    fftlen=\$(( $params.intlen_polcal * 64 ))

    if [ ! -d $baseDir/.deripple_coeffs ]; then
        mkdir $baseDir/.deripple_coeffs
    fi

    args="-f $spectrum"
    args="\$args -l \$fftlen"
    args="\$args -o ${params.label}_polcal_sum_${pol}_f_derippled.npy"
    args="\$args -c $baseDir/.deripple_coeffs"

    python $sourceDir/deripple.py \$args
    """
}

process dedisperse_polcal {
    executor 'slurm'
    cpus 1
    time '5m'
    memory '64 GB'

    input:
    tuple val(pol), path(spectrum) from derippled_spectrum_polcal

    output:
    tuple val(pol), path("${params.label}_polcal_sum_${pol}_f_dedispersed.npy") into dedispersed_spectrum_polcal

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    args="-f $spectrum"
    args="\$args --DM `grep $params.polcal_name $params.cal_file | cut -d " " -f 4`"
    args="\$args --f0 $params.f0_polcal"
    args="\$args --bw 336"
    args="\$args -o ${params.label}_polcal_sum_${pol}_f_dedispersed.npy"

    python3 $sourceDir/dedisperse.py \$args
    """
}

process ifft_polcal {
    executor 'slurm'
    cpus 1
    time '10m'
    memory '32 GB'

    input:
    tuple val(pol), path(spectrum) from dedispersed_spectrum_polcal

    output:
    path("${params.label}_polcal_sum_${pol}_t.npy") into pol_time_series_polcal

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4
    module load scipy/1.4.1-python-3.7.4
    
    args="-f $spectrum"
    args="\$args -o ${params.label}_polcal_sum_${pol}_t.npy"

    python3 $sourceDir/ifft.py \$args
    """
}

process generate_dynspecs_polcal {
    executor 'slurm'
    cpus 1
    time '1h'
    memory '64 GB'

    //publishDir "$sourceDir/output/${params.label}", mode: 'move'

    input:
    path(pol_time_series) from pol_time_series_polcal.collect()

    output:
    path "*dynspec.npy" into dynspecs_polcal

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    args="-x ${params.label}_polcal_sum_x_t.npy"
    args="\$args -y ${params.label}_polcal_sum_y_t.npy"
    args="\$args -o ${params.label}_polcal_sum_!_@.npy"

    python3 $sourceDir/dynspecs.py \$args

    tar -czvhf ${params.label}_polcal_fulltimeres.tar.gz *.npy
    """
}

process polcal {
    publishDir "$baseDir/output/${params.label}", mode: 'copy'

    input:
    path dynspecs_polcal

    output:
    path "*png" into polcal_imgs
    path "polcal_soln.dat" into polcal_soln

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4
    module load astropy/4.0.1-python-3.7.4
    module load scipy/1.6.0-python-3.7.4
    module load matplotlib/3.2.1-python-3.7.4

    args="-i ${params.label}_polcal_sum_i_dynspec.npy"
    args="\$args -q ${params.label}_polcal_sum_q_dynspec.npy"
    args="\$args -u ${params.label}_polcal_sum_u_dynspec.npy"
    args="\$args -v ${params.label}_polcal_sum_v_dynspec.npy"
    args="\$args -p `grep $params.polcal_name $params.cal_file | cut -d " " -f 5`"
    args="\$args -f $params.f0_polcal"
    args="\$args -b 336"
    args="\$args -l ${params.label}"
    args="\$args -r $baseDir/data/vela_parkes/vela_spec_pks_20180920.dat"
    args="\$args -o polcal_soln.dat"
    args="\$args --reduce_df 4"
    args="\$args --plot"
    args="\$args --plotdir ."
    args="\$args --refplotdir $baseDir/data/vela_parkes/"

    python3 $sourceDir/polcal.py \$args
    """
}

process launch_without_polcal {
    output:
    path "polcal_soln.dat" into polcal_soln_skip

    when:
    params.skip_polcal == true

    """
    echo "0" > polcal_soln.dat
    echo "0" >> polcal_soln.dat
    """
}

process beamform_frb {
    executor 'slurm'
    cpus 4
    time '20m'
    memory '32 GB'

    input:
    tuple path(imfile), path(calcfile) from created_calcfiles_frb.mix(got_calcfiles_frb)
    each pol from polarisations_frb
    each antnum from antennas_frb
    path polcal_soln from polcal_soln.mix(polcal_soln_skip)

    output:
    tuple val(pol), path("${params.label}_frb_${antnum}_${pol}_f.npy") into spectra_frb

    """
    module load python/2.7.14 
    module load numpy/1.16.3-python-2.7.14 
    module load scipy/1.0.0-python-2.7.14 
    module load astropy/2.0.3-python-2.7.14 
    module load matplotlib/2.2.2-python-2.7.14 
    module load joblib/0.11

    if [ ! -d delays ]; then
        mkdir delays    # needed by craftcor_tab.py
    fi

    # polarisation calibration solutions
    polcal_delay=`head -1 $polcal_soln`
    polcal_offset=`head -2 $polcal_soln | tail -1`

    args="-i $params.numints_frb"
    args="\$args -n $params.intlen_frb"
    args="\$args --offset $params.offset_frb"
    args="\$args -d $params.data_frb"
    args="\$args --parset $params.fcm"
    args="\$args --calcfile $imfile"
    args="\$args --aips_c $params.bandpass"
    args="\$args --an $antnum"
    args="\$args --pol $pol"
    args="\$args -o ${params.label}_frb_${antnum}_${pol}_f.npy"
    args="\$args --tab"
    args="\$args --polcal_delay=\$polcal_delay"
    args="\$args --polcal_offset=\$polcal_offset"
    if [ ! "${params.hwfile}" = "" ]; then
        args="\$args --hwfile $params.hwfile"
    fi

    echo "python $sourceDir/craftcor_tab.py \$args"
    python $sourceDir/craftcor_tab.py \$args
    """
}

process sum_frb {
    executor 'slurm'
    cpus 1
    time '15m'
    memory '16 GB'

    input:
    tuple val(pol), path(spectra) from spectra_frb.groupTuple()

    output:
    tuple val(pol), path("${params.label}_frb_sum_${pol}_f.npy") into summed_spectrum_frb

    """
    module load gcc/7.3.0 
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    args="--f_dir ."
    args="\$args -f ${params.label}_frb"
    args="\$args -p $pol"
    args="\$args -o ${params.label}_frb_sum_${pol}_f.npy"

    python3 $sourceDir/sum.py \$args
    """
}

process deripple_frb {
    executor 'slurm'
    cpus 1
    time '45m'
    memory '64 GB'

    input:
    tuple val(pol), path(spectrum) from summed_spectrum_frb

    output:
    tuple val(pol), path("${params.label}_frb_sum_${pol}_f_derippled.npy") into derippled_spectrum_frb

    """
    module load python/2.7.14
    module load numpy/1.16.3-python-2.7.14
    module load scipy/1.0.0-python-2.7.14

    fftlen=\$(( $params.intlen_frb * 64 ))

    if [ ! -d $baseDir/.deripple_coeffs ]; then
        mkdir $baseDir/.deripple_coeffs
    fi

    args="-f $spectrum"
    args="\$args -l \$fftlen"
    args="\$args -o ${params.label}_frb_sum_${pol}_f_derippled.npy"
    args="\$args -c $baseDir/.deripple_coeffs"

    python $sourceDir/deripple.py \$args
    """
}

process dedisperse_frb {
    executor 'slurm'
    cpus 1
    time '5m'
    memory '64 GB'

    input:
    tuple val(pol), path(spectrum) from derippled_spectrum_frb

    output:
    tuple val(pol), path("${params.label}_frb_sum_${pol}_f_dedispersed.npy") into dedispersed_spectrum_frb

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    args="-f $spectrum"
    args="\$args --DM $params.DM_frb"
    args="\$args --f0 $params.f0_frb"
    args="\$args --bw 336"
    args="\$args -o ${params.label}_frb_sum_${pol}_f_dedispersed.npy"

    python3 $sourceDir/dedisperse.py \$args
    """
}

process ifft_frb {
    executor 'slurm'
    cpus 1
    time '10m'
    memory '32 GB'

    input:
    tuple val(pol), path(spectrum) from dedispersed_spectrum_frb

    output:
    path("${params.label}_frb_sum_${pol}_t.npy") into pol_time_series_frb

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4
    module load scipy/1.4.1-python-3.7.4
    
    args="-f $spectrum"
    args="\$args -o ${params.label}_frb_sum_${pol}_t.npy"

    python3 $sourceDir/ifft.py \$args
    """
}

process generate_dynspecs_frb {
    executor 'slurm'
    cpus 1
    time '1h'
    memory '64 GB'

    publishDir "$baseDir/output/${params.label}", mode: 'copy'

    input:
    path(pol_time_series) from pol_time_series_frb.collect()

    output:
    path "*npy"

    """
    module load gcc/7.3.0
    module load openmpi/3.0.0
    module load python/3.7.4
    module load numpy/1.18.2-python-3.7.4

    args="-x ${params.label}_frb_sum_x_t.npy"
    args="\$args -y ${params.label}_frb_sum_y_t.npy"
    args="\$args -o ${params.label}_frb_sum_!_@.npy"

    python3 $sourceDir/dynspecs.py \$args

    tar -czvhf ${params.label}_frb_fulltimeres.tar.gz *.npy
    """
}
