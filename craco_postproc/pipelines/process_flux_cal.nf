nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as correlate_fluxcal } from './correlate'
include { determine_flux_cal_solns } from './calibration'

workflow process_flux_cal {
    take:
        label   // val
        target  // val
        data    // val
        polyco  // val
        fcm // val
        ra  // val
        dec // val
        fluxflagfile    // val
        cpasspoly   // val

    main:
        empty_binconfig = create_empty_file("binconfig")
        fits = correlate_fluxcal(
            label, data, fcm, ra, dec, empty_binconfig, polyco, 0, 
            fluxflagfile
        )

        determine_flux_cal_solns(fits, fluxflagfile, target, cpasspoly)
    emit:
        flux_cal_solns = determine_flux_cal_solns.out
}
