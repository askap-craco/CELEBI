nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as correlate_fluxcal } from './correlate'
include { determine_flux_cal_solns } from './calibration'

workflow process_flux_cal {
    take:
        label   // val
        target  // val
        data    // val
        fcm // val
        ra  // val
        dec // val
        cpasspoly   // val
    main:
        empty_binconfig = create_empty_file("binconfig")
        fits = correlate_fluxcal(
            label, data, fcm, ra, dec, empty_binconfig, 0
        )

        println "If you haven't already, you will need to write a flagfile for"
        println "$fits and provide it with --fluxflagfile"

        determine_flux_cal_solns(fits, params.fluxflagfile, target, cpasspoly)
    emit:
        flux_cal_solns = determine_flux_cal_solns.out
}
