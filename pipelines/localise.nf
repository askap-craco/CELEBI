localise_dir = "$baseDir/../localise/"
params.out_dir = "${params.publish_dir}/${params.label}"

process generate_binconfig {
    /*
        Create binconfig files for each correlation mode

        Input
            cand: path
                Candidate to generate binconfigs from

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
    publishDir "${params.out_dir}/binconfigs", mode: "copy"

    input:
        path cand
    output:
        path "craftfrb.finder.binconfig", emit: finder
        path "craftfrb.gate.binconfig", emit: gate
        path "craftfrb.rfi.binconfig", emit: rfi
        path "craftfrb.polyco", emit: polyco
        path "dosubtractions.sh", emit: subtractions
        path "int_time", emit: int_time

    script:
        """
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        # fi
        ml apptainer
        set -a
        set -o allexport
        tmp_file=".TMP_\$BASHPID"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/getGeocentricDelay.py $params.data_frb $cand > \$tmp_file'

        sl2f_cmd=`tail -1 \$tmp_file`
        sl2f_cmd="python3 $localise_dir/\$sl2f_cmd"
        \$sl2f_cmd > sl2f.out
        cat sl2f.out | tail -1 > int_time
        """
    
    stub:
        """
        touch craftfrb.finder.binconfig
        touch craftfrb.gate.binconfig
        touch craftfrb.rfi.binconfig
        touch craftfrb.polyco
        touch dosubtractions.sh
        touch int_time
        """
}

process find_offset {
    /*
        Compare fitted field sources to RACS sources to calculate systematic
        offset in images created from voltages

        Input
            field sources: path
                File containing positions of sources identified in field image
        
        Output
            dat: path
                Files containing RACS source information
	    doff: path
                Details of offsets
            reg: path
                DS9 region file of identified RACS sources
            png: path
                Plots generated while calculating offset for verification and
                troubleshooting
    */
    publishDir "${params.out_dir}/position", mode: "copy"

    input:
        path field_sources
    
    output:
        path "offset0.dat", emit: offset
	path "offsetfit.txt", emit: doffset
        path "*.reg"
        path "*.png"
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            # . $launchDir/../setup_proc
            [ -d  "~/.astropy/cache" ] && rm -r ~/.astropy/cache
        fi
        ml apptainer
        set -a
        set -o allexport
        
        args="-o ${params.label}_RACS.dat"
        args="\$args -a ${params.label}_ASKAP.dat"
        args="\$args -n ${params.label}_names.dat"
        args="\$args -r ${params.label}_RACS_sources.reg"
	args="\$args -j ${params.label}_jmfits.dat"

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/RACS_lookup.py \$args field*jmfit'

        args="--askappos ${params.label}_ASKAP.dat"
        args="\$args --askapnames ${params.label}_names.dat"
	args="\$args --jmfitnames ${params.label}_jmfits.dat"
        args="\$args --fieldfits ${params.out_dir}/finder/${params.label}.fits"
        args="\$args --racs ${params.label}_RACS.dat"
        args="\$args --frbtitletext ${params.label}"

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/src_offsets_rotated.py \$args && python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_rotated_offsets.dat > offsetfit.txt && python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_offsets_unc.dat'
	
	# python3 $localise_dir/weighted_multi_image_fit_updated.py \
        #    askap2racs_rotated_offsets.dat > offsetfit.txt

        #python3 $localise_dir/weighted_multi_image_fit_updated.py \
        #    askap2racs_offsets_unc.dat

        """
    
    stub:
        """
        touch offset0.dat 
	touch offsetfit.txt
        touch stub.reg
        touch stub.png
        """
}

process apply_offset {
    /*
        Apply offset to fitted FRB position

        Input
            offset: path
                Offset as output by weighted_multi_image_fit_updated.py
	    doffset: path
                Detailed offsets
            askap_frb_pos: path
                JMFIT output file of FRB position fit

        Output
            final_position: path
                FRB final position with error as a txt file
	    hpmap: path
                Healpix map in FITS format
    */
    publishDir "${params.out_dir}/position", mode: "copy"
    
    input:
        path offset
	    path doffset
        path askap_frb_pos

    output:
        path "${params.label}_final_position.txt", emit: final_position
	    path "${params.label}_hpmap.FITS", emit: hpmap
    
    script:
        """
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        # fi   
        ml apptainer
        set -a
        set -o allexport
        tmp_file=".TMP_\$BASHPID"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/apply_rotated_offset.py --frbname ${params.label} --frb $askap_frb_pos \
            --offset $offset --doffset $doffset --frbfits ${params.out_dir}/finder/${params.label}.fits \
            --hpfits  ${params.label}_hpmap.FITS > ${params.label}_final_position.txt'        

        """
    
    stub:
        """
        touch ${params.label}_final_position.txt
	    touch ${params.label}_hpmap.FITS
        """
}
