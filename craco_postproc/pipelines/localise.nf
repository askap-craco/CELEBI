#!/usr/bin/env nextflow

/******************************************************************************
 * localisation.nf
 * Author: Danica Scott [danica.scott@postgrad.curtin.edu.au]
 * Last modified: 2021-11-01
 * 
 * This is the pipeline for processing ASKAP CRAFT voltages to obtain FRB
 * localisations
 * 
 * FRB localisation
 * ----------------
 * The FRB localisation pipeline has 19 processes (see localisation_L1.png) 
 * across three broad stages:
 *  1. Correlation
 *      > Correlates three data sets:
 *          - FRB finder: data binned in time to identify the FRB itself
 *          - FRB field: data across the entire duration of the voltages for
 *                       astrometric corrections from field sources
 *          - Calibrator: for calibration of the other two data sets
 *  2. Calibration
 *      > Flag RFI
 *      > Applies calibration solutions determined from the calibrator data set
 *        to the finder and field data sets
 *  3. Source fitting
 *      > Fit the position of the FRB in the finder bin it is detected in
 *      > Identify and fit point sources in the field image, then calculate the
 *        mean offset of the sources by comparing to catalog positions
 *      > Apply the offset to the fitted FRB position to get a final position
 *        and uncertainty region of the FRB
 *
 * Acknowledgements
 * ----------------
 * This file was written by Danica Scott, but is very largely based on the
 * localisation pipeline developed by Adam Deller and Cherie Day (really, it's
 * an automated version of that pipeline).
 *****************************************************************************/

// Parameters that should be set in the user's nextflow.config file
params.user = "dscott"
params.data = "/fred/oz002/users/$params.user/craftpy2/data/"

// Parameters given on the command line
params.calibrate = false
params.calflags = ''
params.frbflags = ''
params.cpasspoly = 3
params.finder_imagesize = 512
params.finder_pixelsize = 1
params.skip_check = false
params.skip_frb = false
params.skip_polcal = false
params.fit = false
params.clean_frb_image = ""
params.frb_pix_coords = ""
params.frb_fit_size = 25
params.clean_field_image = ""
params.field_pix_coords = ""
params.field_fit_size = 25

// Default
params.cal_file = "$baseDir/calibrators.dat"

// Set these in a config file to cut out cards
params.card_list = "1 2 3 4 5 6 7"
params.cards = ['1', '2', '3', '4', '5', '6', '7']
cards = Channel
    .fromList(params.cards)
    .into{cards_cal; cards_polcal; cards_finder; cards_field}

params.fpga_list = "0 1 2 3 4 5"
params.fpgas = ['0', '1', '2', '3', '4', '5']
fpgas = Channel
    .fromList(params.fpgas)
    .into{fpgas_cal; fpgas_polcal; fpgas_finder; 
          fpgas_field}

process check_data {
    /**************************************************************************
     * Sometimes the data download is incomplete. This can be dealt with by
     * either removing affected antennas or by excluding cards. This process 
     * checks if any antennas have incomplete files and stops with a warning to 
     * let you deal with them. This process can be skipped by providing the
     * --skip_check flag.
     *************************************************************************/
    output:
    path "data_checked.out" into data_checked

    when:
    params.skip_check == false

    """
    shortfiles=`ls -l ${params.data_frb}/ak*/* | grep -v hdr | \
                grep -v 29830144 | grep -o ak.._c._f..vcraft`

    echo \$shortfiles

    if [ "\$shortfiles" != "" ]; then
        echo "Some files are short, please check if this is a problem"
        echo "(Skip check with --skip_check)"
        exit 1
    fi

    echo "No short files"
    echo "All good" > data_checked.out
    """
}

// 1.0 - for FRBs
process generate_binconfig {
    /**************************************************************************
     * Create the binconfig files from the detection/snoopy log file. These
     * files define the bins used in the finder branch.
     *************************************************************************/
    output:
    tuple path("craftfrb.finder.binconfig"), path("craftfrb.polyco") into binconfig
    env int_time into int_time_cal, int_time_polcal, int_time_finder, int_time_field

    """
    source $baseDir/setup_difx ${params.user}

    tmp_file=".TMP_\$BASHPID"
    $baseDir/craftpy2/getGeocentricDelay.py ${params.data_frb} ${params.snoopy} > \$tmp_file

    sl2f_cmd=`tail -1 \$tmp_file`
    sl2f_cmd="$baseDir/craftpy2/\$sl2f_cmd"
    \$sl2f_cmd > sl2f.out
    int_time=`cat sl2f.out`
    """
}

process get_startmjd_frb {
    output:
    env startmjd into startmjd_finder, startmjd_field

    """
    module load python/3.8.5

    startmjd=`python3 $baseDir/scripts/get_start_mjd.py $params.data_frb`
    echo "\$startmjd"
    """
}

process get_startmjd_cal {
    output:
    env startmjd into startmjd_cal

    """
    module load python/3.8.5

    startmjd=`python3 $baseDir/scripts/get_start_mjd.py $params.data_cal`
    echo "\$startmjd"
    """
}

process get_startmjd_polcal {
    output:
    env startmjd into startmjd_polcal

    """
    module load python/3.8.5

    startmjd=`python3 $baseDir/scripts/get_start_mjd.py $params.data_polcal`
    echo "\$startmjd"
    """
}

// 2.1
process procTimeStep_cal {
    input:
    env int_time from int_time_cal
    env startmjd from startmjd_cal
    each card from cards_cal
    each f from fpgas_cal

    output:
    path "data_cal/c${card}_f$f" into cal_data
    path "data_cal/c${card}_f$f/*D2D.input" into cal_D2D

    """
    source $baseDir/setup_difx ${params.user}

    # get calibrator RA and dec
    ra_cal=`grep ${params.cal_name} ${params.cal_file} | cut -d " " -f 2`
    dec_cal=`grep ${params.cal_name} ${params.cal_file} | cut -d " " -f 3`

	export CRAFTCATDIR="."

    args="-f ${params.fcm}"
    args="\$args -b 4"
    args="\$args --card $card"
    args="\$args -k"
    args="\$args --slurm"
    args="\$args --name=${params.label}_cal"
    args="\$args -o data_cal"
    args="\$args -t ${params.data_cal}"
    args="\$args --ra \$ra_cal"
    args="\$args -d\$dec_cal"
    args="\$args --freqlabel c${card}_f$f"
    args="\$args --dir=${baseDir}/difx"
    args="\$args --startmjd=\$startmjd"

    $baseDir/craftpy2/processTimeStep.py \$args
    """
}

// 2.1.2
process generate_difx2fits_cal {
    input:
    path cal_data from cal_data.collect()

    output:
    path "data_cal.tar.gz" into all_difx2fits_cal

    """
    mkdir data_cal
    mv c* data_cal/
    cd data_cal
    for c in ${params.card_list}; do
        D2Ds=""
        for f in ${params.fpga_list}; do
            freqlabel="c\${c}_f\${f}"
            D2Dinput=`ls \$freqlabel/*D2D.input`
            D2Ds="\$D2Ds \$D2Dinput"
        done
        difx2fitscmd="difx2fits -v -v -u -B ! \$D2Ds"
        echo "\$difx2fitscmd \"\\\$@\"" | tr ! 0 >> runalldifx2fits
        echo "mv CRAFTFR.0.bin0000.source0000.FITS CRAFT_CARD\$c.FITS" >> runalldifx2fits
    done
    chmod 775 runalldifx2fits
    cd ..
    tar -czvf data_cal.tar.gz data_cal
    """
}

// 2.2
process difx2fits_cal {
    executor 'slurm'
    cpus 1
    memory '4 GB'
    time '10m'

    input:
    path cal_data from all_difx2fits_cal

    output:
    path "CRAFT_CARD*.FITS" into per_card_cal_fits

    """
    source $baseDir/setup_difx ${params.user}

    tar -xzvf $cal_data

    pushd data_cal
    ./runalldifx2fits
    mv CRAFT_CARD*.FITS ..
    popd data_cal
    """
}

// 2.3
process loadfits_cal {
    publishDir "$baseDir/output/${params.label}/${params.cal_name}", mode: "copy"

    input:
    path per_card_cal_fits

    output:
    path "${params.label}_${params.cal_name}.fits" into cal_fits

    """
	echo "source $baseDir/setup_parseltongue ${params.user}"
    source $baseDir/setup_parseltongue ${params.user}

    # get antlist: https://superuser.com/questions/97905/tell-ls-to-print-only-the-base-filename
    antlist=`ls -d ${params.data_cal}/ak* | tr '\\n' '\\0' | xargs -0 -n 1 basename | tr '\\n' ','`

    args="-u 23"
    args="\$args --antlist=\$antlist"
    args="\$args -s 27"
    args="\$args -f ${params.label}_${params.cal_name}.fits"
    args="\$args -o ${params.cal_name}"
    args="\$args CRAFT_CARD?.FITS"

    loadfits.py \$args
    """
}

// 2.1
process procTimeStep_polcal {
    input:
    env int_time from int_time_polcal
    env startmjd from startmjd_polcal
    each card from cards_polcal
    each f from fpgas_polcal

    output:
    path "data_polcal/c${card}_f$f" into polcal_data
    path "data_polcal/c${card}_f$f/*D2D.input" into polcal_D2D

    when:
    params.skip_polcal == false

    """
    source $baseDir/setup_difx ${params.user}

    # get polcalibrator RA and dec
    ra_polcal=`grep ${params.polcal_name} ${params.cal_file} | cut -d " " -f 2`
    dec_polcal=`grep ${params.polcal_name} ${params.cal_file} | cut -d " " -f 3`

	export CRAFTCATDIR="."

    args="-f ${params.fcm}"
    args="\$args -b 4"
    args="\$args --card $card"
    args="\$args -k"
    args="\$args --slurm"
    args="\$args --name=${params.label}_polcal"
    args="\$args -o data_polcal"
    args="\$args -t ${params.data_polcal}"
    args="\$args --ra \$ra_polcal"
    args="\$args -d\$dec_polcal"
    args="\$args --freqlabel c${card}_f$f"
    args="\$args --dir=${baseDir}/difx"
    args="\$args --startmjd=\$startmjd"

    $baseDir/craftpy2/processTimeStep.py \$args
    """
}

// 2.1.2
process generate_difx2fits_polcal {
    input:
    path polcal_data from polcal_data.collect()

    output:
    path "data_polcal.tar.gz" into all_difx2fits_polcal

    """
    mkdir data_polcal
    mv c* data_polcal/
    cd data_polcal
    for c in ${params.card_list}; do
        D2Ds=""
        for f in ${params.fpga_list}; do
            freqlabel="c\${c}_f\${f}"
            D2Dinput=`ls \$freqlabel/*D2D.input`
            D2Ds="\$D2Ds \$D2Dinput"
        done
        difx2fitscmd="difx2fits -v -v -u -B ! \$D2Ds"
        echo "\$difx2fitscmd \"\\\$@\"" | tr ! 0 >> runalldifx2fits
        echo "mv CRAFTFR.0.bin0000.source0000.FITS CRAFT_CARD\$c.FITS" >> runalldifx2fits
    done
    chmod 775 runalldifx2fits
    cd ..
    tar -czvf data_polcal.tar.gz data_polcal
    """
}

// 2.2
process difx2fits_polcal {
    executor 'slurm'
    cpus 1
    memory '4 GB'
    time '10m'

    input:
    path polcal_data from all_difx2fits_polcal

    output:
    path "CRAFT_CARD*.FITS" into per_card_polcal_fits

    """
    source $baseDir/setup_difx ${params.user}

    tar -xzvf $polcal_data

    pushd data_polcal
    ./runalldifx2fits
    mv CRAFT_CARD*.FITS ..
    popd data_polcal
    """
}

// 2.3
process loadfits_polcal {
    publishDir "$baseDir/output/${params.label}/${params.polcal_name}", mode: "copy"

    input:
    path per_card_polcal_fits

    output:
    path "${params.label}_${params.polcal_name}.fits" into polcal_fits

    """
	echo "source $baseDir/setup_parseltongue ${params.user}"
    source $baseDir/setup_parseltongue ${params.user}

    # get antlist: https://superuser.com/questions/97905/tell-ls-to-print-only-the-base-filename
    antlist=`ls -d ${params.data_polcal}/ak* | tr '\\n' '\\0' | xargs -0 -n 1 basename | tr '\\n' ','`

    args="-u 103"
    args="\$args --antlist=\$antlist"
    args="\$args -s 27"
    args="\$args -f ${params.label}_${params.polcal_name}.fits"
    args="\$args -o ${params.polcal_name}"
    args="\$args CRAFT_CARD?.FITS"

    loadfits.py \$args
    """
}

// 3.1
process procTimeStep_finder {
    input:
    tuple path("craftfrb.finder.binconfig"), path("craftfrb.polyco") from binconfig
    env int_time from int_time_finder
    env startmjd from startmjd_finder
    each card from cards_finder
    each f from fpgas_finder

    output:
    path "data_finder/c${card}_f$f" into finder_data
    path "data_finder/c${card}_f$f/*D2D.input" into finder_D2D
    env num_bins into num_bins

    when:
    params.skip_frb == false

    """
    source $baseDir/setup_difx ${params.user}

	export CRAFTCATDIR="."

    args="-f ${params.fcm}"
    args="\$args -b 4"
    args="\$args --card $card"
    args="\$args -k"
    args="\$args --slurm"
    args="\$args --name=${params.label}_finder"
    args="\$args -o data_finder"
    args="\$args -p craftfrb.finder.binconfig"
    args="\$args -t ${params.data_frb}"
    args="\$args -i \$int_time"
    args="\$args --ra ${params.ra_frb}"
    args="\$args -d${params.dec_frb}"
    args="\$args --freqlabel c${card}_f$f"
    args="\$args --dir=${baseDir}/difx"
    args="\$args --startmjd=\$startmjd"

    $baseDir/craftpy2/processTimeStep.py \$args
    
    num_bins=`grep PULSAR craftfrb.finder.binconfig | cut -d ":" -f 2`
    """
}

// 3.1.2
process generate_difx2fits_finder {
    input:
    path finder_data from finder_data.collect()
    val num_bins from num_bins.collect()

    output:
    path "data_finder.tar.gz" into all_difx2fits_finder

    """
    #//num_bins="\${str:num_bins: -3:1}"
    num_bins="${num_bins[0]}"
    echo \$num_bins
    mkdir data_finder
    mv c* data_finder/
    cd data_finder
    for c in ${params.card_list}; do
        D2Ds=""
        for f in ${params.fpga_list}; do
            freqlabel="c\${c}_f\${f}"
            D2Dinput=`ls \$freqlabel/*D2D.input`
            D2Ds="\$D2Ds \$D2Dinput"
        done
        difx2fitscmd="difx2fits -u -B ! \$D2Ds"

        for b in `seq 0 \$((num_bins - 1))`; do
            bin2="\$(printf "%02d" \$b)"
            bin4="\$(printf "%04d" \$b)"
            difx2fitscmd="difx2fits -v -v -u -B \$b \$D2Ds"
            echo "\$difx2fitscmd \"\\\$@\"" >> runalldifx2fits
            echo "mv CRAFTFR.0.bin\$bin4.source0000.FITS CRAFT_CARD\${c}_BIN\${bin2}.FITS" >> runalldifx2fits
        done
    done
    chmod 775 runalldifx2fits
    cd ..
    tar -czvf data_finder.tar.gz data_finder
    """
}

// 3.2
process difx2fits_finder {
    executor 'slurm'
    cpus 1
    memory '4 GB'
    time '10m'

    input:
    path finder_data from all_difx2fits_finder

    output:
    path "CRAFT_CARD*.FITS" into per_card_finder_fits

    """
    source $baseDir/setup_difx ${params.user}

    tar -xzvf $finder_data

    pushd data_finder
    ./runalldifx2fits
    mv CRAFT_CARD*.FITS ..
    popd data_finder
    """
}

// 3.3
process loadfits_finder {
    publishDir "$baseDir/output/${params.label}/finder", mode: "copy"

    input:
    path per_card_finder_fits
    tuple path("craftfrb.finder.binconfig"), path("craftfrb.polyco") from binconfig

    output:
    path "${params.label}_finder_bin*.fits" into finder_fits
    path "loadfits.log" into loadfits_log

    """
    set +e
    source $baseDir/setup_parseltongue ${params.user}

    # get antlist: https://superuser.com/questions/97905/tell-ls-to-print-only-the-base-filename
    antlist=`ls -d ${params.data_cal}/ak* | tr '\\n' '\\0' | xargs -0 -n 1 basename | tr \'\n' ','`

    # get number of bins from binconfig
    num_bins=`grep PULSAR craftfrb.finder.binconfig | cut -d ":" -f 2`

    touch loadfits.log

    for i in `seq 0 \$((num_bins-1))`; do
        bin="\$(printf "%02d" \$i)"

        args="-u 33"
        args="\$args --antlist=\$antlist"
        args="\$args -s 27"
        args="\$args -f ${params.label}_finder_bin\${bin}.fits"
        args="\$args -o FINDB\${bin}"
        args="\$args CRAFT_CARD?_BIN\${bin}.FITS"

        echo "loadfits.py \$args"
        loadfits.py \$args

        # if loadfits didn't work, skip this bin and make a note of it
        if [ ! \$? -eq 0 ]; then
            echo "Bin \$bin failed!"
            echo "loadfits on bin \$bin failed" >> loadfits.log
        fi
    done
    """
}

// 4.1
process procTimeStep_field {
    input:
    env int_time from int_time_field
    env startmjd from startmjd_field
    each card from cards_field
    each f from fpgas_field

    output:
    path "data_field/c${card}_f$f" into field_data
    path "data_field/c${card}_f$f/*D2D.input" into field_D2D

    when:
    params.skip_frb == false

    """
    source $baseDir/setup_difx ${params.user}

    # get a header file for the beam position
    ant_pattern="${params.data_frb}/ak*"
    ants=( \$ant_pattern )
    first_ant=`echo \$ants`
    beam_pattern="\$first_ant/beam*"
    beams=( \$beam_pattern )
    first_beam=`echo \$beams`

    ra_beam_deg=`grep BEAM_RA \$first_beam/*c1_f0*hdr | cut -d " " -f 2`
    dec_beam_deg=`grep BEAM_DEC \$first_beam/*c1_f0*hdr | cut -d " " -f 2`

    echo \$ra_beam_deg
    echo \$dec_beam_deg

    # convert beam RA and dec into hms/dms
    radec_beam=`python $baseDir/get_beam_radec.py \$ra_beam_deg \$dec_beam_deg`
    ra_beam=`echo \$radec_beam | cut -d " " -f 1 | tr h : | tr m : | tr s 0`
    dec_beam=`echo \$radec_beam | cut -d " " -f 2 | tr d : | tr m : | tr s 0`


	export CRAFTCATDIR="."

    args="-f ${params.fcm}"
    args="\$args -b 4"
    args="\$args --card $card"
    args="\$args -k"
    args="\$args --slurm"
    args="\$args --name=${params.label}_field"
    args="\$args -o data_field"
    args="\$args -i \$int_time"
    args="\$args -t ${params.data_frb}"
    args="\$args --ra \$ra_beam"
    args="\$args -d\$dec_beam"
    args="\$args --freqlabel c${card}_f$f"
    args="\$args --dir=${baseDir}/difx"
    args="\$args --startmjd=\$startmjd"

    $baseDir/craftpy2/processTimeStep.py \$args
    """
}

// 4.1.2
process generate_difx2fits_field {
    input:
    path field_data from field_data.collect()

    output:
    path "data_field.tar.gz" into all_difx2fits_field

    """
    mkdir data_field
    mv c* data_field/
    cd data_field
    for c in ${params.card_list}; do
        D2Ds=""
        for f in ${params.fpga_list}; do
            freqlabel="c\${c}_f\${f}"
            D2Dinput=`ls \$freqlabel/*D2D.input`
            D2Ds="\$D2Ds \$D2Dinput"
        done
        difx2fitscmd="difx2fits -v -v -u -B ! \$D2Ds"
        echo "\$difx2fitscmd \"\\\$@\"" | tr ! 0 >> runalldifx2fits
        echo "mv CRAFTFR.0.bin0000.source0000.FITS CRAFT_CARD\$c.FITS" >> runalldifx2fits
    done
    chmod 775 runalldifx2fits
    cd ..

    tar -czvf data_field.tar.gz data_field
    """
}

// 4.2
process difx2fits_field {
    executor 'slurm'
    cpus 1
    memory '4 GB'
    time '10m'

    input:
    path field_data from all_difx2fits_field

    output:
    path "CRAFT_CARD*.FITS" into per_card_field_fits

    """
    source $baseDir/setup_difx ${params.user}

    tar -xzvf $field_data

    pushd data_field
    ./runalldifx2fits
    mv CRAFT_CARD*.FITS ..
    popd data_field
    """
}

// 4.3
process loadfits_field {
    publishDir "$baseDir/output/${params.label}/field", mode: "copy"

    input:
    path per_card_field_fits

    output:
    path "${params.target}_field.fits" into field_fits

    """
    source $baseDir/setup_parseltongue ${params.user}

    # get antlist: https://superuser.com/questions/97905/tell-ls-to-print-only-the-base-filename
    antlist=`ls -d ${params.data_frb}/ak* | tr '\\n' '\\0' | xargs -0 -n 1 basename | tr '\\n' ','`

    args="-u 43"
    args="\$args --antlist=\$antlist"
    args="\$args -s 27"
    args="\$args -f ${params.label}_field.fits"
    args="\$args -o ${params.target}"
    args="\$args CRAFT_CARD?.FITS"

    loadfits.py \$args
    """
}

// 5.1
process calibrate_cal {
    publishDir "$baseDir/output/${params.label}/calibrated/$params.cal_name", mode: "copy"

    input:
    path cal_fits

    output:
    tuple path("calibration_noxpol_${params.target}.tar.gz"), path("${params.label}_${params.cal_name}_calibrated_uv.fits"), path("${params.label}_${params.cal_name}_calibrated_uv.ms") into cal_solns_finder, cal_solns_field, cal_solns_ips, cal_solns_polcal
    path "*noxpol*" into bandpass

    when:
    params.calibrate == true

    """
    source $baseDir/setup_parseltongue ${params.user}

    args="--calibrateonly"
    args="\$args -c $cal_fits"
    args="\$args --uvsrt"
    args="\$args -u 51"
    args="\$args --src=${params.target}"
    args="\$args --cpasspoly=${params.cpasspoly}"
    args="\$args -f 15"
    args="\$args --flagfile=${params.calflags}"

    calibrateFRB.py \$args
    """
}

process calibrate_polcal {
    publishDir "$baseDir/output/${params.label}/calibrated/$params.polcal_name", mode: "copy"

    input:
    path polcal_fits
    tuple path(cal_solns), path(cal_fits), path(cal_ms) from cal_solns_polcal

    output:
    path "${params.label}_polcal_calibrated.tar.gz" into calibrated_polcal_fits

    """
    source $baseDir/setup_parseltongue ${params.user}

    tar -xzvf $cal_solns

    args="--targetonly"
    args="\$args -t $polcal_fits"
    args="\$args -r 3"
    args="\$args --cpasspoly=${params.cpasspoly}"
    args="\$args -i"
    args="\$args --dirtymfs"
    args="\$args -a 16"
    args="\$args -u 500"
    args="\$args --skipplot"
    args="\$args --tarflagfile=${params.polcalflags}"
    args="\$args --src=${params.target}"

    calibrateFRB.py \$args
    
    tar -czvf ${params.label}_polcal_calibrated.tar.gz *calibrated*
    """
}

// 5.2.2
process calibrate_finder {
    publishDir "$baseDir/output/${params.label}/calibrated/finder", mode: "copy"

    input:
    each finder_fits
    tuple path(cal_solns), path(cal_fits), path(cal_ms) from cal_solns_finder

    output:
    path "${params.label}_finder_calibrated_bin*.tar.gz" into calibrated_finder_fits

    when:
    params.skip_frb == false

    """
    source $baseDir/setup_parseltongue ${params.user}

    tar -xzvf $cal_solns

    bin=`echo $finder_fits | cut -d '_' -f 3 | cut -d 'n' -f 2 | cut -d '.' -f 1`
    mkdir bin\$bin

    args="--targetonly"
    args="\$args -t $finder_fits"
    args="\$args -r 3"
    args="\$args --cpasspoly=${params.cpasspoly}"
    args="\$args -u 52\$bin"
    args="\$args -i"
    args="\$args --dirtymfs"
    args="\$args -a 16"
    args="\$args --imagesize=${params.finder_imagesize}"
    args="\$args --pixelsize=${params.finder_pixelsize}"
    args="\$args --src=${params.target}"
    args="\$args --skipplot"
    args="\$args --tarflagfile=${params.frbflags}"

    calibrateFRB.py \$args

    mv TARGET* bin\$bin
    mv *calibrated* bin\$bin

    tar -czvf ${params.label}_finder_calibrated_bin\${bin}.tar.gz bin\$bin
    """
}

// 5.3
process calibrate_field {
    publishDir "$baseDir/output/${params.label}/calibrated/field", mode: "copy"

    input:
    path field_fits
    tuple path(cal_solns), path(cal_fits), path(cal_ms) from cal_solns_field

    output:
    path "${params.label}_field_calibrated.tar.gz" into calibrated_field_fits

    when:
    params.skip_frb == false

    """
    source $baseDir/setup_parseltongue ${params.user}

    tar -xzvf $cal_solns

    args="--targetonly"
    args="\$args -t $field_fits"
    args="\$args -r 3"
    args="\$args --cpasspoly=${params.cpasspoly}"
    args="\$args -u 53"
    args="\$args --skipplot"
    args="\$args --tarflagfile=${params.frbflags}"
    args="\$args --src=${params.target}"

    calibrateFRB.py \$args
    
    tar -czvf ${params.label}_field_calibrated.tar.gz *calibrated*
    """
}

process fit_FRB_position {
    output:
    path "${params.label}_frb*" into frb_fit

    when:
    params.fit == true
    params.clean_frb_image != ""
    params.frb_pix_coords != ""

    """
    . setup_parseltongue $params.user

    args="-f $params.clean_frb_image"
    args="\$args -c $params.frb_pix_coords"
    args="\$args -o ${params.label}_frb"
    args="\$args -s $params.frb_fit_size"
    args="\$args -a ${params.label}_frb_askap.dat"

    $baseDir/scripts/do_jmfit.py \$args
    """
}

process fit_field_sources {
    output:
    path "${params.label}_field*" into field_fit

    when:
    params.fit == true
    params.clean_field_image != ""
    params.field_pix_coords != ""

    """
    . setup_parseltongue $params.user

    args="-f $params.clean_field_image"
    args="\$args -c $params.field_pix_coords"
    args="\$args -o ${params.label}_field"
    args="\$args -s $params.field_fit_size"
    args="\$args -a ${params.label}_field_askap.dat"

    $baseDir/scripts/do_jmfit.py \$args
    """
}

process apply_offset {
    input:
    path frb_fit
    path field_fit

    """

    """
}
