nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as correlate_fluxcal } from './correlate'
include { determine_flux_cal_solns } from './calibration'
include { generate_binconfig } from './localise'

workflow process_flux_cal {
    take:
        label   // val
        target  // val
        data    // val
        data_frb // val
        snoopy // val
        fcm // val
        ra  // val
        dec // val
        fluxflagfile    // val
        cpasspoly   // val

    main:
        if( params.nocorrelate ) {
            fits = Channel.fromPath("${params.publish_dir}/${params.label}/loadfits/fluxcal/*fits")
        }
        else {
            binconfig = generate_binconfig(data_frb, snoopy)
            empty_binconfig = create_empty_file("binconfig")
            fits = correlate_fluxcal(
                label, data, fcm, ra, dec, empty_binconfig, binconfig.polyco, empty_binconfig, 
                fluxflagfile, "fluxcal"
            )
        }

        if( params.calibrate ) {
            flux_cal_solns = determine_flux_cal_solns(fits, fluxflagfile, target, cpasspoly).solns
        }
        else {
            flux_cal_solns = ""
        }
    emit:
        flux_cal_solns
}
