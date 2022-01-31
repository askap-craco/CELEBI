nextflow.enable.dsl=2

include { create_empty_file } from 'utils'
include { correlate } from 'correlate'

workflow process_flux_cal {
    take:
        val label
        val data
        val fcm
        val ra
        val dec
    main:
        fits = correlate(
            label, data, fcm, ra, dec, create_empty_file("binconfig"), 0
        )
    emit:
        flux_cal_solns = determine_flux_cal_solns.out
}
