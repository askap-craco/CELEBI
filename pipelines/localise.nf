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

    container "file://$params.container"

    label 'python'

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
        source /opt/setup_proc_container 
        set -xu

        tmp_file=".TMP_\$BASHPID"
        python3 $localise_dir/getGeocentricDelay.py $params.data_frb $cand > \$tmp_file

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

    container "file://$params.container"

    label 'python'

    input:
        path field_sources
    
    output:
        path "offset0.dat", emit: offset
	path "offsetfit.txt", emit: doffset
        path "*.reg"
        path "*.png"
    
    script:
        """
        source /opt/setup_proc_container
        set -xu

        hostname >> hostname.txt
        python3 $localise_dir/RACS_lookup.py \
               -o ${params.label}_RACS.dat \
               -a ${params.label}_ASKAP.dat \
               -n ${params.label}_names.dat \
               -r ${params.label}_RACS_sources.reg \
	           -j ${params.label}_jmfits.dat \
               field*jmfit

        python3 $localise_dir/src_offsets_rotated.py \
                --askappos ${params.label}_ASKAP.dat \
                --askapnames ${params.label}_names.dat \
	            --jmfitnames ${params.label}_jmfits.dat \
                --fieldfits ${params.out_dir}/finder/${params.label}.fits \
                --racs ${params.label}_RACS.dat \
                --frbtitletext ${params.label}

        python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_rotated_offsets.dat > offsetfit.txt 
        python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_offsets_unc.dat
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
    
    container "file://$params.container"

    label 'python'

    input:
        path offset
	    path doffset
        path askap_frb_pos

    output:
        path "${params.label}_final_position.txt", emit: final_position
    
    script:
        """
        source /opt/setup_proc_container
        set -xu
        
        tmp_file=".TMP_\$BASHPID"
        python3 $localise_dir/apply_rotated_offset.py \
                --frbname ${params.label} \
                --frb $askap_frb_pos \
                --offset $offset \
                --doffset $doffset \
                --frbfits ${params.out_dir}/finder/${params.label}.fits  \
                > ${params.label}_final_position.txt
        """

    
    stub:
        """
        touch ${params.label}_final_position.txt
	    touch ${params.label}_hpmap.FITS
        """
}
