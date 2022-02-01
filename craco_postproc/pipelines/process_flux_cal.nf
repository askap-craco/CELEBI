nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate } from './correlate'

process determine_flux_cal_solns {
    input:
        path cal_fits
        path flagfile
        val target
        val cpasspoly

    output:
        path "calibration_noxpol_${target}.tar.gz"

    script:
        """
        args="--calibrateonly"
        args="\$args -c $cal_fits"
        args="\$args --uvsrt"
        args="\$args -u 51"
        args="\$args --src=$target"
        args="\$args --cpasspoly=$cpasspoly"
        args="\$args -f 15"
        args="\$args --flagfile=$flagfile"

        calibrateFRB.py \$args
        """
}

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
