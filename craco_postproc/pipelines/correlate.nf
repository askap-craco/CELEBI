// Cards and FPGAs to be processed. Override these in a config file to cut out
// data
params.cards = ["1", "2", "3", "4", "5", "6", "7"]
cards = Channel.fromList(params.cards)
params.fpgas = ["0", "1", "2", "3", "4", "5"]
fpgas = Channel.fromList(params.fpgas)

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
    output:
    env startmjd

    """
    startmjd=`python $baseDir/scripts/get_start_mjd.py $params.data
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
    input:
    env startmjd
    each card from cards
    each fpga from fpgas

    output:
    path "c${card}_f${fpga}" into correlated_data
    path "c${card}_f${fpga}/*D2D.input" into D2D

    """
    export CRAFTCATDIR="."  # necessary?

    args="-f $params.fcm"
    args="\$args -b 4"
    args="\$args -k"
    args="\$args --name=$params.label"
    args="\$args -o ."
    args="\$args -t $params.data"
    args="\$args --ra $params.ra"
    args="\$args -d$params.dec"
    args="\$args --card $card"
    args="\$args --freqlabel c${card}_f${fpga}"
    args="\$args --dir=$baseDir/difx"
    args="\$args --startmjd=\$startmjd"

    if [ "$params.binconfig" != "" ]; then
        args="\$args -p $params.binconfig"
    fi
    if [ "$params.inttime" != "" ]; then
        args="\$args -i $params.inttime"
    fi

    $baseDir/craco-postproc/localise/processTimeStep.py \$args
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
    path correlated_data from correlated_data.collect()

    output:
    path "*.FITS" into per_card_fits

    """
    for c in `seq 1 7`; do
        D2Ds=""
        dirs=`find c\$c* 2> /dev/null`
        for d in \$dirs; do
           D2Dinput=`ls \$dir/*D2D.input`
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

    *******************************************************************/
    input:
    path per_card_fits

    output:
    path "${params.label}.fits" into fits

    """
    antlist=`ls -d ${params.data}/ak* | tr '\\n' '\\0' | xargs -0 -n 1 basename | tr '\\n' ','`

    args="-u 1"
    args="\$args --antlist=\$antlist"
    args="\$args -s 27"
    args="\$args -f ${params.label}.fits"
    args="\$args -o ${params.label}"
    args="\$args CRAFT_CARD?.FITS"

    loadfits.py \$args
    """
}
