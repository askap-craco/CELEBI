localise_dir = "$projectDir/../localise"
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
        if [ "$params.uselocalracs" == "true" ]; then
            args="\$args --localracssourcepath=${params.localracssourcepath}"
        fi

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && hostname >> hostname.txt && python3 $localise_dir/RACS_lookup.py \$args field*jmfit'

        args="--askappos ${params.label}_ASKAP.dat"
        args="\$args --askapnames ${params.label}_names.dat"
        args="\$args --jmfitnames ${params.label}_jmfits.dat"
        args="\$args --fieldfits ${params.out_dir}/finder/${params.label}.fits"
        args="\$args --racs ${params.label}_RACS.dat"
        args="\$args --frbtitletext ${params.label}"

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/src_offsets_rotated.py \$args && python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_rotated_offsets.dat > offsetfit.txt && python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_offsets_unc.dat'

        mkdir noexclusions
        mv ${params.label}_ASKAP.dat noexclusions/
        mv ${params.label}_RACS.dat noexclusions/
        mv ${params.label}_names.dat noexclusions/
        mv ${params.label}_RACS_sources.dat noexclusions/
        mv ${params.label}_jmfits.dat noexclusions/
        mv ${params.label}_field_offsets_from_racs.png noexclusions/
        mv askap2racs_offsets_unc.dat noexclusions/
        mv askap2racs_rotated_offsets.dat noexclusions/
        mv offsetfit.txt noexclusions/
        mv offset0.dat noexclusions/
        mv err_vs_offset_0.pdf noexclusions/

        # Re-run with exclusions
        args="-o ${params.label}_RACS.dat"
        args="\$args -a ${params.label}_ASKAP.dat"
        args="\$args -n ${params.label}_names.dat"
        args="\$args -r ${params.label}_RACS_sources.reg"
        args="\$args -j ${params.label}_jmfits.dat"
        if [ "$params.uselocalracs" == "true" ]; then
            args="\$args --localracsgausspath=${params.localracsgausspath}"
            args="\$args --localracssourcepath=${params.localracssourcepath}"
        fi

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/RACS_lookup.py \$args field*jmfit'

        args="--askappos ${params.label}_ASKAP.dat"
        args="\$args --askapnames ${params.label}_names.dat"
        args="\$args --jmfitnames ${params.label}_jmfits.dat"
        args="\$args --fieldfits ${params.out_dir}/finder/${params.label}.fits"
        args="\$args --racs ${params.label}_RACS.dat"
        args="\$args --frbtitletext ${params.label}"

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/src_offsets_rotated.py \$args && python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_rotated_offsets.dat > offsetfit.txt && python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_offsets_unc.dat'
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
    
    script:
        """
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        # fi   
        ml apptainer
        set -a
        set -o allexport
        tmp_file=".TMP_\$BASHPID"

        # Get the FRB galactic latitude
        if [ "$params.uselocalracs" == "true" ]; then
            apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/frb_galactic_coordinates.py --frbra=${params.ra_frb} --frbdec=${params.dec_frb} --planelatcut ${params.racs_plane_latcut} --onplanera ${params.racs_rasystematics_onplane} --offplanera ${params.racs_rasystematics_offplane} --onplanedec ${params.racs_decsystematics_onplane} --offplanedec ${params.racs_decsystematics_offplane} > frameuncertainties.txt'
        fi 

        args="--frbname ${params.label}"
        args="\$args --frb $askap_frb_pos"
        args="\$args --offset $offset"
        args="\$args --doffset $doffset"
        args="\$args --frbfits ${params.out_dir}/finder/${params.label}.fits"
        if [ "$params.uselocalracs" == "true" ]; then
            rasys=`awk '{print \$1}' frameuncertainties.txt`
            decsys=`awk '{print \$2}' frameuncertainties.txt`
            args="\$args --framerauncertainty=\$rasys"
            args="\$args --framedecuncertainty=\$decsys"
        fi

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/apply_rotated_offset.py \$args > ${params.label}_final_position.txt'
        """

    
    stub:
        """
        touch ${params.label}_final_position.txt
	    touch ${params.label}_hpmap.FITS
        """
}
