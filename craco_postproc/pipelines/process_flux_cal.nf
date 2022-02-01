nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate } from './correlate'
include { determine_flux_cal_solns } from './calibration'

workflow process_flux_cal {
    take:
        val label
        val target
        val data
        val fcm
        val ra
        val dec
        val cpasspoly
    main:
        empty_binconfig = create_empty_file("binconfig")
        fits = correlate(
            label, data, fcm, ra, dec, empty_binconfig, 0
        )

        if( params.flagfile )
            determine_flux_cal_solns(fits, params.flagfile, target, cpasspoly)
        else
            println "Please write flagfile for $fits!"
    emit:
        flux_cal_solns = determine_flux_cal_solns.out
}
