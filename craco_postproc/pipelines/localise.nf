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
        python3 $localise_dir/getGeocentricDelay.py $data $snoopy > \$tmp_file

        sl2f_cmd=`tail -1 \$tmp_file`
        sl2f_cmd="python3 $localise_dir/\$sl2f_cmd"
        \$sl2f_cmd > sl2f.out
        int_time=`cat sl2f.out`
        """
}

process apply_offset {
    publishDir "${params.publish_dir}/${params.label}/position", mode: "copy"

    input:
        path field_sources
        path askap_frb_pos
    
    output:
        path "${params.label}_final_position.txt" 
        path "*.dat" 
        path "*.reg"
        path "*.png"
    
    script:
        """
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

        python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_offsets_unc.dat

        python3 $localise_dir/apply_offset.py --frb $askap_frb_pos --offset offset0.dat > ${params.label}_final_position.txt
        """
}
