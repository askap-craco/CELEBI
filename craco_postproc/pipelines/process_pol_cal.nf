nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as correlate_polcal } from './correlate'
include { beamform as beamform_polcal } from './beamform'
include { apply_flux_cal_solns_polcal; determine_pol_cal_solns } from './calibration'
include { generate_binconfig } from './localise'

workflow process_pol_cal {
    take:
        label   // val
        target  // val
        data    // val
        data_frb
        snoopy
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
        binconfig = generate_binconfig(data_frb, snoopy)
        empty_file = create_empty_file("file")
        fits = correlate_polcal(
            label, data, fcm, ra0, dec0, empty_file, binconfig.polyco, 0, polflagfile, "polcal"
        )

        pos = apply_flux_cal_solns_polcal(
            fits, flux_cal_solns, polflagfile, target, cpasspoly
        ).jmfit

        htr_data = beamform_polcal(
            label, data, fcm, pos, flux_cal_solns, empty_file,
            num_ints, int_len, offset, dm, centre_freq
        )
        determine_pol_cal_solns(htr_data)           

    emit:
        pol_cal_solns = determine_pol_cal_solns.out.pol_cal_solns
}
