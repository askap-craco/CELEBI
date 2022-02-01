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
