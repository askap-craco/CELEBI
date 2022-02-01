nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate } from './correlate'
include { beamform } from './beamform'
include { apply_flux_cal_solns, determine_pol_cal_solns } from './calibration'

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

workflow process_pol_cal {
    take:
        label   // val
        target  // val
        data    // val
        fcm // val
        ra0 // val
        dec0    // val
        cpasspoly   // val
        flux_cal_solns  // path
        num_ints    // val
        int_len // val
        offset  // val
        dm  // val
        centre_freq // val

    main:
        empty_binconfig = create_empty_file("binconfig")
        fits = correlate(
            label, data, fcm, ra0, dec0, empty_binconfig, 0
        )

        if( params.flagfile )
            image = apply_flux_cal_solns(
                fits, flux_cal_solns, params.flagfile, target, cpasspoly
            )
            pos = localise(image)
            empty_pol_cal_solns = create_empty_file("polcal.dat")
            htr_data = beamform(
                label, data, fcm, pos, flux_cal_solns, empty_pol_cal_solns, 
                num_ints, int_len, offset, dm, centre_freq
            )
            determine_pol_cal_solns(htr_data)
        else
            println "Please write flagfile for $fits!"            

    emit:
        pol_cal_solns = determine_pol_cal_solns.out
}
