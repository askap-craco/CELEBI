#!/usr/bin/env nextflow
nextflow.enable.dsl=2   // Enable DSL2

// Cards and FPGAs to be processed. Override these in a config file to cut out
// data. The lowest card-fpga pair is used as a reference correlation.
params.cards = ["1", "2", "3", "4", "5", "6", "7"]
cards = Channel.fromList(params.cards)
params.fpgas = ["0", "1", "2", "3", "4", "5"]
fpgas = Channel.fromList(params.fpgas)
card_fpgas = cards.combine(fpgas)
    .filter{ !(it[0] == params.cards.min() & it[1] == params.fpgas.min()) }
ref_card_fpga = cards.min().combine(fpgas.min())

localise_dir = "$projectDir/../localise"

params.uppersideband = false
params.out_dir = "${params.publish_dir}/${params.label}"

process get_start_mjd {
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
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        # fi
        ml apptainer
        set -a
        set -o allexport
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/get_start_mjd.py $data' 
 
        """
    
    stub:
        """
        echo "0"
        """
}

process do_ref_correlation {
    /*
        Correlate the data using DiFX

        Input            
            label: val
                FRB name and context of process instance as a string (no 
                spaces)
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
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
        val ra
        val dec
        path binconfig, stageAs: "craftfrb.binconfig"
        path polyco, stageAs: "craftfrb.polyco"
        path inttime
        val startmjd
        tuple val(card), val(fpga)
        path fcm

    output:
        path "c${card}_f${fpga}", emit: cx_fy
        path "c${card}_f${fpga}/*D2D.input"

    script:
        """
        export CRAFTCATDIR="."  # necessary?
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        # fi

        # create .bat0
        ml apptainer
        set -a
        set -o allexport
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && bat0.pl `find $data/*/*/*vcraft | head -1`'
        # bat0.pl `find $data/*/*/*vcraft | head -1`

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
        args="\$args --dir=$projectDir/../difx"
        args="\$args --startmjd=$startmjd"

        # High-band FRBs need --uppersideband
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
        int_time=`cat $inttime`
        if [ "\$int_time" != "" ]; then
            args="\$args -i \$int_time"
        fi
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && echo \$args && python3 $localise_dir/processTimeStep.py \$args'
        """
    
    stub:
        """
        mkdir c${card}_f${fpga}
        touch c${card}_f${fpga}/stubD2D.input
        """
}
process do_correlation {
    /*
        Correlate the data using DiFX

        Input            
            label: val
                FRB name and context of process instance as a string (no 
                spaces)
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
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
            refcorr: path
                Optional reference correlation for fillDiFX. If not needed,
                pass an empty file
            fcm: path
                fcm file to use
        
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
        val ra
        val dec
        path binconfig, stageAs: "craftfrb.binconfig"
        path polyco, stageAs: "craftfrb.polyco"
        path inttime
        val startmjd
        tuple path(ref_corr), val(card), val(fpga)
        path fcm

    output:
        path "c${card}_f${fpga}", emit: cx_fy
        path "c${card}_f${fpga}/*D2D.input"

    script:
        """
        export CRAFTCATDIR="."  # necessary?
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        # fi
        ml apptainer
        set -a
        set -o allexport
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && bat0.pl `find $data/*/*/*vcraft | head -1`'
        # create .bat0
        # bat0.pl `find $data/*/*/*vcraft | head -1`

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
        args="\$args --dir=$projectDir/../difx"
        args="\$args --startmjd=$startmjd"

        # High-band FRBs need --uppersideband
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
        int_time=`cat $inttime`
        if [ "\$int_time" != "" ]; then
            args="\$args -i \$int_time"
        fi

        # Provide reference correlation
        args="\$args --ref=$ref_corr"

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && echo \$args && python3 $localise_dir/processTimeStep.py \$args'

        """
    
    stub:
        """
        mkdir c${card}_f${fpga}
        touch c${card}_f${fpga}/stubD2D.input
        """
}

process difx_to_fits {
    /*
        Convert correlated data from the DiFX format into FITS.

        Input
            label: val
                FRB name and context of process instance as a string (no 
                spaces)
            correlated_data: path
                All the c*_f* directories produced by do_correlation
            mode: val
                If mode == "finder", loop over finder bins to do conversion,
                otherwise do it without looping
        Output:
            fits: path
                A single FITS file containing data across all cards processed
    */
    publishDir "${params.out_dir}/loadfits/${mode}", mode: "copy"

    input:
        val label
        path correlated_data
        path polyco
        val mode

    output:
        path "${label}*.fits", emit: fits
        path "finderbin04.fits", emit: centre, optional: true

    script:
        """
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        # fi
        ml apptainer
        set -a
        set -o allexport
        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"

        for c in `seq 1 7`; do
            D2Ds=""
            if find c\$c*; then
                dirs=`find c\$c* 2> /dev/null`
                for d in \$dirs; do
                    D2Dinput=`ls \$d/*D2D.input`
                    D2Ds="\$D2Ds \$D2Dinput"
                done

                if [ "$mode" == "finder" ]; then
                    for b in `seq 0 7`; do
                        bin2="\$(printf "%02d" \$b)"
                        bin4="\$(printf "%04d" \$b)"
                        difx2fitscmd="difx2fits -v -v -u -B \$b \$D2Ds"
                        echo "\$difx2fitscmd \"\\\$@\"" | tr ! 0 >> runalldifx2fits
                        echo "mv CRAFTFRB.0.bin\${bin4}.source0000.FITS \
                            CRAFT_CARD\${c}_BIN\${bin2}.FITS" >> runalldifx2fits
                    done
                else
                    difx2fitscmd="difx2fits -v -v -u -B ! \$D2Ds"
                    echo "\$difx2fitscmd \"\\\$@\"" | tr ! 0 >> runalldifx2fits
                    echo "mv CRAFTFRB.0.bin0000.source0000.FITS \
                        CRAFT_CARD\$c.FITS" >> runalldifx2fits
                fi
            fi
        done
        chmod 775 runalldifx2fits
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && ./runalldifx2fits'

        antlist=""
        for i in `seq -w 1 36`; do
            antlist="\${antlist}ak\$i,"
        done
        echo \$antlist

        label=$label
        label=\${label:0:12}    # Truncate label to fit in AIPS

        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"

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

            # if [ "$params.ozstar" == "true" ]; then
            #    echo ". $launchDir/../setup_loadfits
            # fi
            echo "loadfits.py \$args" >> doloadfits
            chmod 775 doloadfits
            apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ./doloadfits'
        else
            for i in `seq 0 7`; do
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
                # if [ "$params.ozstar" == "true" ]; then
                #    echo ". $launchDir/../setup_parseltongue3" | tr ! 0 >> doloadfits
                # fi    
                # echo ". $launchDir/../setup_parseltongue3" | tr ! 0 >> doloadfits
                echo "loadfits.py \$args" >> doloadfits
                chmod 775 doloadfits
                apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ./doloadfits'
                rm -f doloadfits
            done
        fi
        rm -rf aips_dir
        """
    
    stub:
        """
        if [ "$mode" == "finder" ]; then
            for i in `seq 0 7`; do
                touch finderbin0\${i}.fits
            done
        else
            touch ${label}_stub.fits
        fi
        """
}

process subtract_rfi {
    /*
        Subtract RFI visibilities from target visibilities to remove RFI from
        the data without zapping channels that may contain very important
        signal!

        Input
            target_fits: path
                A single target bin's visibility FITS file
            rfi_fits: path
                RFI-only visibility FITS file
            subtractions: path
                file containing subtractions commands with correctly calculated
                scale argument. If an empty file, assume a scale of 1
                
        Output
            fits: path
                Visbility FITS file for a target bin with RFI subtracted
    */
    maxForks 1

    input:
        each path(target_fits)
        path rfi_fits
        path subtractions
    
    output:
        path "*.fits"
    
    script:
        """
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_parseltongue3
        # fi

        # subtractions not empty: finder mode
        ml apptainer
        set -a
        set -o allexport
        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"
        if [ `wc -c $subtractions | awk '{print \$1}'` != 0 ]; then
            fits="$target_fits"
            bin=\${fits:9:2}
            sleep \$bin     # stagger starts of parallel processes
            scale=\$(grep finderbin00.fits dosubtractions.sh | cut -d' ' -f4)

            apptainer exec $params.container bash -c 'source /opt/setup_proc_container && uvsubScaled.py $target_fits *_rfi.fits \$scale norfifbin\${bin}.fits'
        else
            apptainer exec $params.container bash -c 'source /opt/setup_proc_container && uvsubScaled.py $target_fits $rfi_fits 1 norfi_$target_fits'
        fi
        rm -rf aips_dir
        """
    
    stub:
        """
        fits="$target_fits"
        bin=\${fits:9:2}
        touch norfifbin\${bin}.fits
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
            ra: val
                Right ascension to correlate around, as hh:mm:ss
            dec: val
                Declination to correlate around, as dd:mm:ss
            binconfig: path
                File specifying how to bin the correlated data
            polyco: path
                TODO: describe polyco
            inttime: path
                File containing integration time in seconds
            mode: val
                String describing specific mode (i.e. finder, field, rfi)
            fcm: path
                fcm to use
        
        Emit
            fits: path
                FITS visibility file correlated from voltages
    */
    take:
        label
        data
        ra
        dec
        binconfig
        polyco
        inttime
        mode
        fcm

    main:
        startmjd = get_start_mjd(data)
        println(data)
        // reference correlation
        ref_correlation = do_ref_correlation(
            label, data, ra, dec, binconfig, polyco, inttime, startmjd, 
            ref_card_fpga, fcm
        ).cx_fy

        // card_fpgas kicks off an instance of do_correlation for
        // every unique card-fpga pair, which are then collated with .collect()
        correlated_data = do_correlation(
            label, data, ra, dec, 
            binconfig.first(), 
            polyco.first(), 
            inttime.first(), 
            startmjd, 
            ref_correlation.combine(card_fpgas),
            fcm
        ).cx_fy

        all_correlations = ref_correlation.concat(correlated_data).collect()

        (fits, centre) = difx_to_fits(
            label, all_correlations, polyco, mode
        )
    
    emit:
        fits = fits
        centre = centre
}
