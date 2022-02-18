#!/usr/bin/env nextflow
nextflow.enable.dsl=2   // Enable DSL2

// Cards and FPGAs to be processed. Override these in a config file to cut out
// data
params.cards = ["1", "2", "3", "4", "5", "6", "7"]
cards = Channel.fromList(params.cards)
params.fpgas = ["0", "1", "2", "3", "4", "5"]
fpgas = Channel.fromList(params.fpgas)

localise_dir = "$baseDir/../localise/"

process get_startmjd {
    /*******************************************************************
    Find the start time of the data by finding the earliest listed 
    startmjd in the data's headers. This is to prevent edge cases where 
    files start on either side of an integer second boundary, which 
    results in mismatched data down the line.

    Output:
    startmjd - environment variable
        The earliest start time found in the data headers
    *******************************************************************/
    input:
    val data
    
    output:
    stdout

    """
    python $localise_dir/get_start_mjd.py $data
    """
}

process process_time_step {
    /*******************************************************************
    Correlate the data using DiFX

    Input:
    startmjd - environment variable
        The earliest start time found in the data headers
    card - value
        Specific card to be processed by this instance
    fpga - value
        Specific FPGA to be processed by this instance
    
    Output:
    correlated_data - path
        A directory "c{card}_f{fpga}" containing the correlated data
    D2D - path
        A file within the correlated_data directory that only exists if 
        the correlation has completed. This is used as an output to 
        ensure that Nextflow correctly determines if this process has 
        completed or not, since the correlated_data directory will exist
        even if the process has started without completing.
    *******************************************************************/
    cpus 15
    
    input:
    val label
    val data
    val fcm
    val ra
    val dec
    path binconfig, stageAs: "craftfrb.binconfig"
    path polyco, stageAs: "craftfrb.polyco"
    val inttime
    val startmjd
    tuple val(card), val(fpga)

    output:
    path "c${card}_f${fpga}", emit: cx_fy
    path "c${card}_f${fpga}/*D2D.input"

    """
    export CRAFTCATDIR="."  # necessary?

    args="-f $fcm"
    args="\$args -b 4"
    args="\$args -k"
    args="\$args --name=$label"
    args="\$args -o ."
    args="\$args -t $data"
    args="\$args --ra $ra"
    args="\$args -d$dec"
    args="\$args --card $card"
    freqlabel="c${card}_f${fpga}"
    args="\$args --freqlabel \$freqlabel"
    args="\$args --dir=$baseDir/../difx"
    args="\$args --startmjd=$startmjd"
    args="\$args --ts=5"

    mkdir \$freqlabel
    cp craftfrb.polyco \$freqlabel

    # Only use binconfig if it's not empty
    if [ `wc -c craftfrb.binconfig | awk '{print \$1}'` != 0 ]; then
        args="\$args -p craftfrb.binconfig"
    fi

    # Only include inttime if non-zero
    if [ "$inttime" != "0" ]; then
        args="\$args -i $inttime"
    fi

    $localise_dir/processTimeStep.py \$args
    """
}

process difx2fits {
    /*******************************************************************
    Convert the correlated data from the DiFX format into FITS.

    Input:
    correlated_data - collection of paths
        All the c*_f* directories produced by process_time_step
    
    Output:
    per_card_fits - collection of paths
        A .FITS visibilities file for each card processed
    *******************************************************************/
    input:
    path correlated_data

    output:
    path "*.FITS"

    """
    for c in `seq 1 7`; do
        D2Ds=""
        dirs=`find c\$c* 2> /dev/null`
        for d in \$dirs; do
           D2Dinput=`ls \$d/*D2D.input`
           D2Ds="\$D2Ds \$D2Dinput"
        done
        difx2fitscmd="difx2fits -v -v -u -B ! \$D2Ds"
        echo "\$difx2fitscmd \"\\\$@\"" | tr ! 0 >> runalldifx2fits
        echo "mv CRAFTFR.0.bin0000.source0000.FITS CRAFT_CARD\$c.FITS" >> runalldifx2fits
    done
    chmod 775 runalldifx2fits
    ./runalldifx2fits
    """

}

process loadfits {
    /*******************************************************************
    Combine the per-card FITS files into a single FITS file

    Input:
    per_card_fits - collection of paths
        A FITS file for each card processed
    
    Output:
    fits - path
        A single FITS file containing data across all cards processed
    *******************************************************************/
    input:
    val data
    val label
    path per_card_fits
    val flagfile

    output:
    path "${label}.fits"

    script:
        """
        antlist=`ls -d $data/ak* | tr '\\n' '\\0' | xargs -0 -n 1 basename | tr '\\n' ','`

        label=$label
        label=\${label:0:12}    # Truncate label to fit in AIPS

        args="-u 1"
        args="\$args --antlist=\$antlist"
        args="\$args -s 27"
        args="\$args -f ${label}.fits"
        args="\$args -o \$label"
        args="\$args CRAFT_CARD?.FITS"

        loadfits.py \$args

        if [ "$flagfile" == "" ]; then
            echo "You now need to write the flagfile for ${label}.fits!"
            exit 2
        """
}

process subtract_rfi {
    input:
        path target_fits
        path rfi_fits
    
    output:
        path "rfi_subtracted.fits"
    
    script:
        """
        #TODO: is the scale factor here a constant, or should it be calculated?
        uvsubScaled.py $target_fits $rfi_fits 0.159380579 rfi_subtracted.fits
        """
}

workflow correlate {
    take:
        label   // val
        data    // val
        fcm // val
        ra  // val
        dec // val
        binconfig   // path
        polyco  // path
        inttime // val
        flagfile    // val

    main:
        startmjd = get_startmjd(data)

        // cards.combine(fpgas) kicks off an instance of process_time_step for
        // every unique card-fpga pair, which are then collated with .collect()
        correlated_data = process_time_step(
            label, data, fcm, ra, dec, binconfig, polyco, 0, startmjd, 
            cards.combine(fpgas)
        )
        per_card_fits = difx2fits(correlated_data.cx_fy.collect())

        loadfits(data, label, per_card_fits, flagfile)
    
    emit:
        fits = loadfits.out
}
