#!/usr/bin/env nextflow
nextflow.enable.dsl=2   // Enable DSL2

// Cards and FPGAs to be processed. Override these in a config file to cut out
// data
params.cards = ["1", "2", "3", "4", "5", "6", "7"]
cards = Channel.fromList(params.cards)
params.fpgas = ["0", "1", "2", "3", "4", "5"]
fpgas = Channel.fromList(params.fpgas)

localise_dir = "$baseDir/../localise/"

params.uppersideband = false

process get_startmjd {
    /*
        Find the start time of the data by finding the earliest listed 
        startmjd in the data's headers. This is to prevent edge cases where 
        files start on either side of an integer second boundary, which 
        results in mismatched data down the line.

        Input
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)

        Output
            stdout
                The earliest start time found in the data headers in MJD
    */
    input:
        val data
    
    output:
        stdout

    script:
        """    
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi
        python3 $localise_dir/get_start_mjd.py $data
        """
}

process create_bat0 {
    /*
        Create .bat0 file from vcraft data.
        TODO: describe what the .bat0 file contains

        Input
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
        
        Output
            bat0: path
                .bat0 file
    */
    input:
        val data

    output:
        path ".bat0", emit: bat0

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi
        bat0.pl `find $data/*/*/*vcraft | head -1`
        """
}

process process_time_step {
    /*
        Correlate the data using DiFX

        Input            
            label: val
                FRB name and context of process instance as a string (no 
                spaces)
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
            fcm: val
                Absolute path to fcm (hardware delays) file
            ra: val
                Right ascension to correlate around, as hh:mm:ss
            dec: val
                Declination to correlate around, as dd:mm:ss
            binconfig: path
                File specifying how to bin the correlated data
            polyco: path
                TODO: describe polyco
            inttime: val
                Integration time in seconds
            startmjd: environment variable
                The earliest start time found in the data headers
            card, fpga: tuple(val, val)
                Specific card-fpga pair to be correlated by this instance
            bat0: path
                TODO: describe bat0
        
        Output
            correlated_data: path
                A directory "c{card}_f{fpga}" containing the correlated data
            D2D: path
                A file within the correlated_data directory that only exists if 
                the correlation has completed. This is used as an output to 
                ensure that Nextflow correctly determines if this process has 
                completed or not, since the correlated_data directory will 
                exist even if the process has ended (e.g. by being killed) 
                without completing.
    */
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
        path bat0

    output:
        path "c${card}_f${fpga}", emit: cx_fy
        path "c${card}_f${fpga}/*D2D.input"

    script:
        """
        export CRAFTCATDIR="."  # necessary?
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi

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

        # High-band FRBs need --upersideband
        if [ "$params.uppersideband" == "true" ]; then
            args="\$args --uppersideband"
        fi

        # if running on ozstar, use the slurm queue, otherwise run locally 
        # across 16 cpus
        if [ "$params.ozstar" == "true" ]; then
            args="\$args --slurm"
        else
            args="\$args --ts 16"
        fi

        if [ -d \$freqlabel ]; then
            rm -r \$freqlabel
        fi
        mkdir \$freqlabel
        cp craftfrb.polyco \$freqlabel
        cp .bat0 \$freqlabel

        # Only use binconfig if it's not empty
        if [ `wc -c craftfrb.binconfig | awk '{print \$1}'` != 0 ]; then
            args="\$args -p craftfrb.binconfig"
        fi

        # Only include inttime if non-zero
        if [ "$inttime" != "0" ]; then
            args="\$args -i $inttime"
        fi

        python3 $localise_dir/processTimeStep.py \$args
        """
}

process difx2fits {
    /*
        Convert correlated data from the DiFX format into FITS.

        Input
            correlated_data: path
                All the c*_f* directories produced by process_time_step
            polyco: path
                TODO: describe polyco 
            mode: val
                If mode == "finder", loop over finder bins to do conversion,
                otherwise do it without looping
        
        Output
            per_card_fits: path
                A .FITS visibilities file for each card processed
    */
    input:
        path correlated_data
        path polyco, stageAs: "craftfrb.polyco"
        val mode

    output:
        path "*.FITS"

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi
        for c in `seq 1 7`; do
            D2Ds=""
            if find c\$c*; then
                dirs=`find c\$c* 2> /dev/null`
                for d in \$dirs; do
                    D2Dinput=`ls \$d/*D2D.input`
                    D2Ds="\$D2Ds \$D2Dinput"
                done

                if [ "$mode" == "finder" ]; then
                    for b in `seq 0 20`; do
                        bin2="\$(printf "%02d" \$b)"
                        bin4="\$(printf "%04d" \$b)"
                        difx2fitscmd="difx2fits -v -v -u -B \$b \$D2Ds"
                        echo "\$difx2fitscmd \"\\\$@\"" | tr ! 0 >> runalldifx2fits
                        echo "mv CRAFTFR.0.bin\${bin4}.source0000.FITS \
                            CRAFT_CARD\${c}_BIN\${bin2}.FITS" >> runalldifx2fits
                    done
                else
                    difx2fitscmd="difx2fits -v -v -u -B ! \$D2Ds"
                    echo "\$difx2fitscmd \"\\\$@\"" | tr ! 0 >> runalldifx2fits
                    echo "mv CRAFTFR.0.bin0000.source0000.FITS \
                        CRAFT_CARD\$c.FITS" >> runalldifx2fits
                fi
            fi
        done
        chmod 775 runalldifx2fits
        ./runalldifx2fits
        """

}

process loadfits {
    /*
        Combine the per-card FITS files into a single FITS file

        Input            
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
            label: val
                FRB name and context of process instance as a string (no 
                spaces)
                TODO: data and label are in the reverse order to elsewhere
            per_card_fits: path
                A FITS file for each card processed
            mode: val
                If mode == "finder", loop over finder bins, otherwise do it 
                without looping
        
        Output:
            fits: path
                A single FITS file containing data across all cards processed
    */
    publishDir "${params.publish_dir}/${params.label}/loadfits/${mode}", mode: "copy"

    input:
        val data
        val label
        path per_card_fits
        val mode

    output:
        path "${label}*.fits"

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi
        antlist=`ls -d $data/ak* | tr '\\n' '\\0' | xargs -0 -n 1 basename | tr '\\n' ','`

        label=$label
        label=\${label:0:12}    # Truncate label to fit in AIPS

        aipsid="\$((RANDOM%8192))"
        if [ "$mode" != "finder" ]; then
            args="-u \$aipsid"
            #args="-u \${BASHPID: -4}"   # get randomly-generated user id
            args="\$args --antlist=\$antlist"
            args="\$args -s 27"
            args="\$args -f ${label}.fits"
            args="\$args -o \$label"
            args="\$args CRAFT_CARD?.FITS"

            echo "loadfits.py \$args"
            #loadfits.py \$args

            if [ "$params.ozstar" == "true" ]; then
                echo ". $launchDir/../setup_proc" | tr ! 0 >> doloadfits
            fi
            echo "loadfits.py \$args" >> doloadfits
            chmod 775 doloadfits
            ./doloadfits
        else
            for i in `seq 0 20`; do
                bin="\$(printf "%02d" \$i)"
                args="-u \$aipsid"
                #args="-u \${BASHPID: -4}"   # get randomly-generated user id
                args="\$args --antlist=\$antlist"
                args="\$args -s 27"
                args="\$args -f ${label}bin\${bin}.fits"
                args="\$args -o \${label}bin\${bin}"
                args="\$args CRAFT_CARD?_BIN\${bin}.FITS"

                echo "loadfits.py \$args"
                #loadfits.py \$args
                if [ "$params.ozstar" == "true" ]; then
                    echo ". $launchDir/../setup_proc" | tr ! 0 >> doloadfits
                fi    
                echo ". $launchDir/../setup_proc" | tr ! 0 >> doloadfits
                echo "loadfits.py \$args" >> doloadfits
                chmod 775 doloadfits
                ./doloadfits
                rm -f doloadfits
            done
        fi
        """
}

process subtract_rfi {
    /*
        Subtract RFI visibilities from finder visibilities to remove RFI from
        the data without zapping channels that may contain very important
        signal!

        Input
            finder_fits: path
                A single finder bin's visibility FITS file
            rfi_fits: path
                RFI-only visibility FITS file
            subtractions: path
                file containing subtractions commands with correctly calculated
                scale argument
                
        Output
            fits: path
                Visbility FITS file for a finder bin with RFI subtracted
    */
    maxForks 1

    input:
        each path(finder_fits)
        path rfi_fits
        path subtractions
    
    output:
        path "*.fits"
    
    script:
        """
        fits="$finder_fits"
        bin=\${fits:9:2}
        sleep \$bin     # stagger starts of parallel processes
        scale=\$(grep finderbin00.fits dosubtractions.sh | cut -d' ' -f4)

        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi

        uvsubScaled.py $finder_fits *_rfi.fits \$scale norfifbin\${bin}.fits
        """
}

workflow correlate {
    /*
        Workflow to correlate vcraft voltages into visibilities

        Take
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
            fcm: val
                Absolute path to fcm (hardware delays) file
            ra: val
                Right ascension to correlate around, as hh:mm:ss
            dec: val
                Declination to correlate around, as dd:mm:ss
            binconfig: path
                File specifying how to bin the correlated data
            polyco: path
                TODO: describe polyco
            inttime: val
                Integration time in seconds
            mode: val
                String describing specific mode (i.e. finder, field, rfi)
        
        Emit
            fits: path
                FITS visibility file correlated from voltages
    */
    take:
        label
        data
        fcm
        ra
        dec
        binconfig
        polyco
        inttime
        mode

    main:
        startmjd = get_startmjd(data)
        bat0 = create_bat0(data)

        // cards.combine(fpgas) kicks off an instance of process_time_step for
        // every unique card-fpga pair, which are then collated with .collect()
        correlated_data = process_time_step(
            label, data, fcm, ra, dec, binconfig, polyco, inttime, startmjd, 
            cards.combine(fpgas), bat0
        )
        per_card_fits = difx2fits(correlated_data.cx_fy.collect(), polyco, mode)

        loadfits(data, label, per_card_fits, mode)
    
    emit:
        fits = loadfits.out
}
