localise_dir = "$baseDir/../localise/"

process generate_binconfig {
    /*
        Create binconfig files for each correlation mode

        Input
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
        
        Output
            TODO: describe the modes here
            finder: path
                Finder mode binconfig
            gate: path
                Gate mode binconfig
            rfi: path
                RFI mode binconfig
            polyco: path
                TODO: describe polyco
            subtractions: path
                File containing subtractions commands with correctly calculated
                scale argument
            int_time: env
                Integration time in seconds
    */
    input:
        val data
        path snoopy

    output:
        path "craftfrb.finder.binconfig", emit: finder
        path "craftfrb.gate.binconfig", emit: gate
        path "craftfrb.rfi.binconfig", emit: rfi
        path "craftfrb.polyco", emit: polyco
        path "dosubtractions.sh", emit: subtractions
        env int_time, emit: int_time

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi
        tmp_file=".TMP_\$BASHPID"
        python3 $localise_dir/getGeocentricDelay.py $data $snoopy > \$tmp_file

        sl2f_cmd=`tail -1 \$tmp_file`
        sl2f_cmd="python3 $localise_dir/\$sl2f_cmd"
        \$sl2f_cmd > sl2f.out
        int_time=`cat sl2f.out | tail -1`
        """
}

process apply_offset {
    /*
        Compare fitted field sources to RACS sources to calculate systematic
        offset in images created from voltages, then apply that offset to the
        fitted FRB position to obtain the final localisation

        Input
            field sources: path
                File containing positions of sources identified in field image
            askap_frb_pos: path
                JMFIT output file of FRB position fit
        
        Output
            final_position: path
                Text file containing final FRB position with errors
            dat: path
                Files containing RACS source information
            reg: path
                DS9 region file of identified RACS sources
            png: path
                Plots generated while calculating offset for verification and
                troubleshooting
    */
    publishDir "${params.publish_dir}/${params.label}/position", mode: "copy"

    input:
        path field_sources
        path askap_frb_pos
    
    output:
        path "${params.label}_final_position.txt", emit: final_position
        path "*.dat" 
        path "*.reg"
        path "*.png"
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
            [ -d  "~/.astropy/cache" ] && rm -r ~/.astropy/cache
        fi        
        args="-o ${params.label}_RACS.dat"
        args="\$args -a ${params.label}_ASKAP.dat"
        args="\$args -n ${params.label}_names.dat"
        args="\$args -r ${params.label}_RACS_sources.reg"

        python3 $localise_dir/RACS_lookup.py \$args field*jmfit

        args="--askappos ${params.label}_ASKAP.dat"
        args="\$args --askapnames ${params.label}_names.dat"
        args="\$args --racs ${params.label}_RACS.dat"
        args="\$args --frbtitletext ${params.label}"

        python3 $localise_dir/src_offsets.py \$args

        python3 $localise_dir/weighted_multi_image_fit_updated.py \
            askap2racs_offsets_unc.dat

        python3 $localise_dir/apply_offset.py --frb $askap_frb_pos \
            --offset offset0.dat > ${params.label}_final_position.txt
        """
}
