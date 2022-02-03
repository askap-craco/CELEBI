nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as correlate_polcal } from './correlate'
include { beamform as beamform_polcal } from './beamform'
include { apply_flux_cal_solns; determine_pol_cal_solns } from './calibration'
include { localise as localise_polcal } from './localise'

workflow process_pol_cal {
    take:
        label   // val
        target  // val
        data    // val
        fcm // val
        ra0 // val
        dec0    // val
        polflagfile // val
        cpasspoly   // val
        flux_cal_solns  // path
        num_ints    // val
        int_len // val
        offset  // val
        dm  // val
        centre_freq // val

    main:
        empty_binconfig = create_empty_file("binconfig")
        fits = correlate_polcal(
            label, data, fcm, ra0, dec0, empty_binconfig, 0, polflagfile
        )

        image = apply_flux_cal_solns(
            fits, flux_cal_solns, polflagfile, target, cpasspoly
        )
        pos = localise_polcal(image)
        empty_pol_cal_solns = create_empty_file("polcal.dat")
        htr_data = beamform_polcal(
            label, data, fcm, pos, flux_cal_solns, empty_pol_cal_solns,
            num_ints, int_len, offset, dm, centre_freq
        )
        determine_pol_cal_solns(htr_data)           

    emit:
        pol_cal_solns = determine_pol_cal_solns.out
}
