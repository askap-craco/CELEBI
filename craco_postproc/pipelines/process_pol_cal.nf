nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate } from './correlate'
include { beamform } from './beamform'
include { apply_flux_cal_solns, determine_pol_cal_solns } from './calibration'
include { localise } from './localise'

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
