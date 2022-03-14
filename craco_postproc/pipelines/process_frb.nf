nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as correlate_finder; correlate as correlate_rfi;
    correlate as correlate_field; subtract_rfi as subtract_rfi_finder;
    subtract_rfi as subtract_rfi_field } from './correlate'
include { apply_flux_cal_solns_finder;
    apply_flux_cal_solns_field } from './calibration'
include { localise; apply_offset; generate_binconfig } from './localise'
include { beamform as beamform_frb } from './beamform'

workflow process_frb {
    take:
        label   // val
        data    // val
        snoopy
        fcm // val
        ra0 // val
        dec0    // val
        fieldflagfile
        flux_cal_solns  // path
        // pol_cal_solns // path
        cpasspoly   // val
        num_ints    // val
        int_len // val
        offset  // val
        dm  // val
        centre_freq // val

    main:
        binconfig = generate_binconfig(data, snoopy)
        binconfig_finder = binconfig.finder
        binconfig_rfi = binconfig.rfi
        subtractions = binconfig.subtractions
        polyco = binconfig.polyco
        int_time = binconfig.int_time

        finder_fits = correlate_finder(
            "finder", data, fcm, ra0, dec0, binconfig_finder, polyco, int_time, "N/A", "finder"
        )
        rfi_fits = correlate_rfi(
            "${label}_rfi", data, fcm, ra0, dec0, binconfig_rfi, polyco, int_time, "N/A", "rfi"
        )
        empty_file = create_empty_file("file")
        field_fits = correlate_field(
            "${label}_field", data, fcm, ra0, dec0, empty_file, polyco, 0, "N/A", "rfi"
        )

        no_rfi_finder_fits = subtract_rfi_finder(finder_fits, rfi_fits, subtractions, "finder")

        finder_image = apply_flux_cal_solns_finder(
            no_rfi_finder_fits, flux_cal_solns, label, cpasspoly
        ).peak_image

        apply_flux_cal_solns_field(
            field_fits, flux_cal_solns, fieldflagfile, label, cpasspoly, finder_image
        )

        // apply_offset(field_image, askap_frb_pos)

        // beamform_frb(
        //     label, data, fcm, askap_frb_pos, flux_cal_solns, pol_cal_solns,
        //     num_ints, int_len, offset, dm, centre_freq
        // )

    // emit:
        // true_pos = apply_offset.out
        // frb_htr_data = beamform_frb.out
}
