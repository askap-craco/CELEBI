nextflow.enable.dsl=2   // Enable DSL2
nextflow.enable.strict=true // be less generous

flagging_dir = "${projectDir}/../flagging/"
utils_dir    = "${projectDir}/../utils/"
beamform_dir = "${projectDir}/../beamform/"
localise_dir = "${projectDir}/../localise/"

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

    label 'celebi'

    input:
        path cand
    output:
        path "craftfrb.finder.binconfig", emit: finder
        path "craftfrb.gate.binconfig", emit: gate
        path "craftfrb.rfi.binconfig", emit: rfi
        path "craftfrb.polyco", emit: polyco
        path "dosubtractions.sh", emit: subtractions
        path "int_time", emit: int_time
        path "geo_delay.txt", emit: geo_delay

    script:
        """
        source /opt/setup_proc_container 
        set -xu
        
        tmp_file=".TMP_\$BASHPID"
        python3 $localise_dir/getGeocentricDelay.py $params.data_frb $cand $params.numfinderbins $params.searchms > \$tmp_file

        python3 $localise_dir/\$(tail -1 \$tmp_file) > sl2f.out

        tail -1 sl2f.out  > int_time
        """
    
    stub:
        """
        touch craftfrb.finder.binconfig
        touch craftfrb.gate.binconfig
        touch craftfrb.rfi.binconfig
        touch craftfrb.polyco
        touch dosubtractions.sh
        touch int_time
        touch geo_delay.txt
        """
}

process find_offset {
    /*
        Compare fitted field sources to RACS sources to calculate systematic
        offset in images created from voltages

        Input
            field fits: path
                FITS image for the field
            field sources: path
                File containing positions of sources identified in field image
            exlabel: val
                Exclusive label        
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

    label 'celebi'

    input:
        path fld_fits
        path field_sources
        val exlabel
    
    output:
        path "offset0_*.dat", emit: offset
	    path "offsetfit_*.txt", emit: doffset
        path "*.reg"
        path "*.png"
    
    script:
        """
        source /opt/setup_proc_container
        set -xu

        hostname >> hostname.txt
        
        racsvlass=' '
        if [ "$params.uselocalcatalog" == "true" ]; then
            if [ "$params.referencecatalog" == "RACS" ]; then 
                racsvlass="--localracssourcepath=${params.localracssourcepath}"
            elif [ "$params.referencecatalog" == "VLASS" ]; then
                racsvlass="--localvlasspath=${params.localvlasspath}"
            else
                echo "Not sure what to do with reference catalog ${params.referencecatalog}"
            fi 
        fi

        python3 $localise_dir/getmatchradius.py ${fld_fits} ${params.matchradius}
        
        radius_arcsec=\$(cat "match_radius_arcsec.txt")
        echo "Matching radius (arcsec) = "\${radius_arcsec}

        python3 $localise_dir/RACS_lookup.py \
               -o ${params.label}_RACS.dat \
               -a ${params.label}_ASKAP.dat \
               -n ${params.label}_names.dat \
               -r ${params.label}_RACS_sources.reg \
	           -j ${params.label}_jmfits.dat \
               --matchrad=\${radius_arcsec} \
	           --referencecatalog=${params.referencecatalog} \
	           \$racsvlass \
               field*jmfit

        python3 $localise_dir/src_offsets_rotated.py \
                --askappos ${params.label}_ASKAP.dat \
                --askapnames ${params.label}_names.dat \
	            --jmfitnames ${params.label}_jmfits.dat \
                --fieldfits ${fld_fits} \
                --racs ${params.label}_RACS.dat \
                --frbtitletext ${params.label}

        python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_rotated_offsets.dat > offsetfit.txt 
        python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_offsets_unc.dat
        
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
        
        racsvlass=' '
        if [ "$params.uselocalcatalog" == "true" ]; then
            if [ "$params.referencecatalog" == "RACS" ]; then 
                racsvlass="--localracsgausspath=${params.localracsgausspath} --localracssourcepath=${params.localracssourcepath}"
            elif [ "$params.referencecatalog" == "VLASS" ]; then
                racsvlass="--localvlasspath=${params.localvlasspath}"
            else
                echo "Not sure what to do with reference catalog ${params.referencecatalog}"
            fi 
        fi
        
        python3 $localise_dir/RACS_lookup.py \
               -o ${params.label}_RACS.dat \
               -a ${params.label}_ASKAP.dat \
               -n ${params.label}_names.dat \
               -r ${params.label}_RACS_sources.reg \
	           -j ${params.label}_jmfits.dat \
               --matchrad=\${radius_arcsec} \
	           --referencecatalog=${params.referencecatalog} \
	           \$racsvlass \
               field*jmfit
        
        python3 $localise_dir/src_offsets_rotated.py \
                --askappos ${params.label}_ASKAP.dat \
                --askapnames ${params.label}_names.dat \
	            --jmfitnames ${params.label}_jmfits.dat \
                --fieldfits ${fld_fits} \
                --racs ${params.label}_RACS.dat \
                --frbtitletext ${params.label}
        
        python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_rotated_offsets.dat > offsetfit_${exlabel}.txt 
        python3 $localise_dir/weighted_multi_image_fit_updated.py askap2racs_offsets_unc.dat

        mv offset0.dat offset0_${exlabel}.dat

        """
    
    stub:
        """
        touch offset0.dat 
        touch offsetfit.txt
        touch stub.reg
        touch stub.png
        touch offset0_${exlabel}.dat
        touch offsetfit_${exlabel}.txt
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
            exlabel: val
                Excluseive label  
        Output
            final_position: path
                FRB final position with error as a txt file
            hpmap: path
                Healpix map in FITS format
    */
    publishDir "${params.out_dir}/position", mode: "copy"

    label 'celebi'
    // label 'conda'

    input:
        path offset
	    path doffset
        path askap_frb_pos
        val exlabel

    output:
        path "*_final_position.txt", emit: final_position
    
    script:
        """
        source /opt/setup_proc_container
        set -xu

        tmp_file=".TMP_\$BASHPID"
        
        # Get the FRB galactic latitude
        if [ "$params.uselocalcatalog" == "true" ]; then
            python3 $localise_dir/frb_galactic_coordinates.py --frbra=${params.ra_frb} --frbdec=${params.dec_frb} --planelatcut=${params.catalog_plane_latcut} \
                    --onplanera=${params.catalog_rasystematics_onplane} --offplanera=${params.catalog_rasystematics_offplane} \
                    --onplanedec=${params.catalog_decsystematics_onplane} --offplanedec=${params.catalog_decsystematics_offplane} > frameuncertainties.txt
        fi
        
        radecsys=' '
        if [ "$params.uselocalcatalog" == "true" ]; then
            rasys=`awk '{print \$1}' frameuncertainties.txt`
            decsys=`awk '{print \$2}' frameuncertainties.txt`
            radecsys="--framerauncertainty=\$rasys --framedecuncertainty=\$decsys"
        fi
        
        python3 $localise_dir/apply_rotated_offset.py --frbname=${params.label} --frb=${askap_frb_pos} --offset=$offset --doffset=$doffset \
                --frbfits=${params.out_dir}/finder/${params.label}.fits \$radecsys > ${params.label}_${exlabel}_final_position.txt

        """
    
    stub:
        """
        touch ${params.label}_final_position.txt
	    touch ${params.label}_hpmap.FITS
        """
}

process find_frb_beam_position {
    /*
        Find the FRB position relative to the beam

        Input
            askap_frb_pos: path
                JMFIT output file of FRB position fit
            fieldfits: path
                FITS image of the field  
        Output
            Beam info: path
                Text file with all relevant info
            Position plot: path
                Plot of FRB position w.r.t. the primary beam
    */
    publishDir "${params.out_dir}/position", mode: "copy"

    label 'celebi'
    // label 'conda'

    input:
        path frbposfile
	    path fieldfits

    output:
        path "*.txt", emit: beaminfo
        path "*.png", emit: beampos
    
    script:
        """
        source /opt/setup_proc_container
        set -xu

        vcfile=\$(find ${params.data_frb} -type f -name "ak*_c1_f1.vcraft.hdr" | head -n 1)
        
        python3 $localise_dir/find_frb_beam_info.py --posfile $frbposfile \
            --refvcraftfile \$vcfile \
            --fieldfits $fieldfits \
            --cenfreqmhz ${params.centre_freq_frb} \
            --bwmhz ${params.bw} \
            --diametre ${params.askapdishdia} \
            --plotpng frb_beam_pos > frb_beam_info.txt

        """
    
    stub:
        """
        touch frb_beam_info.txt
	    touch frb_beam_pos.png
        """
}
