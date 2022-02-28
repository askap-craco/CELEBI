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
    input:
        path field_image
        tuple val(ra_askap), val(dec_askap)
    
    output:
        tuple val(ra_true), val(dec_true)
    
    exec:
        ra_true = 0
        dec_true = 0
}
