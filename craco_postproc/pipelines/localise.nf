localise_dir = "$baseDir/../localise/"

process generate_binconfig {
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
        tmp_file=".TMP_\$BASHPID"
        $localise_dir/getGeocentricDelay.py $data $snoopy > \$tmp_file

        sl2f_cmd=`tail -1 \$tmp_file`
        sl2f_cmd="$localise_dir/\$sl2f_cmd"
        \$sl2f_cmd > sl2f.out
        int_time=`cat sl2f.out`
        """
}

process localise {
    input:
        path image

    output:
        tuple val(ra), val(dec)

    script:
        """
        localise $image > pos.dat   # temp
        """
    
    // Mixing script block with exec block might not work, but trying anyway
    exec:
        // format of pos.dat:
        //  RA (hms)
        //  uRA (arcsec)
        //  Dec (dms)
        //  uDec (arcsec)
        reader = file('pos.dat').newReader()
        ra = reader.readLine()
        ura = reader.readLine()
        dec = reader.readLine()
        udec = reader.readLine()
}

process apply_offset {
    publishDir "${params.publish_dir}/${params.label}", mode: "copy"

    input:
        path field_sources
        path askap_frb_pos
    
    output:
        path "${params.label}_final_position.txt" 
        path "*.dat" 
        path "*.reg" 
    
    script:
        """
        args="-o ${params.label}_RACS.dat"
        args="\$args -a ${params.label}_ASKAP.dat"
        args="\$args -n ${params.label}_names.dat"
        args="\$args -r ${params.label}_RACS_sources.reg"

        $localise_dir/RACS_lookup.py \$args field*jmfit

        args="--askappos ${params.label}_ASKAP.dat"
        args="\$args --askapnames ${params.label}_names.dat"
        args="\$args --racs ${params.label}_RACS.dat"
        args="\$args --frbtitletext ${params.label}"

        $localise_dir/src_offsets.py \$args

        $localise_dir/weighted_multi_image_fit_updated.py askap2racs_offsets_unc.dat

        $localise_dir/apply_offset.py --frb $askap_frb_pos --offset offset0.dat > ${params.label}_final_position.txt
        """
}
