nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_finder; correlate as corr_rfi;
    correlate as corr_field; subtract_rfi as subtract_rfi_finder } from './correlate'
include { apply_flux_cal_solns_finder as cal_finder;
    apply_flux_cal_solns_field as cal_field; get_peak } from './calibration'
include { apply_offset; generate_binconfig } from './localise'
include { beamform as bform_frb } from './beamform'

params.fieldimage = ""
params.flagfinder = false

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
        pol_cal_solns // path
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
        if ( new File("${params.publish_dir}/${params.label}/loadfits/finder/finderbin20.fits").exists() ) {
            finder_fits = Channel.fromPath("${params.publish_dir}/${params.label}/loadfits/finder/*fits")
        }
        else {
            finder_fits = corr_finder(
                "finder", data, fcm, ra0, dec0, binconfig_finder, polyco, int_time, "N/A", "finder"
            )
        }

        if ( new File("${params.publish_dir}/${params.label}/loadfits/rfi/${params.label}_rfi.fits").exists() ) {
            rfi_fits = Channel.fromPath("${params.publish_dir}/${params.label}/loadfits/rfi/${params.label}_rfi.fits")
        }
        else {
            rfi_fits = corr_rfi(
                "${label}_rfi", data, fcm, ra0, dec0, binconfig_rfi, polyco, int_time, "N/A", "rfi"
            )
        }

        if ( (params.fieldimage != "") or new File("${params.publish_dir}/${params.label}/loadfits/field/${params.label}_field.fits").exists() ) {
            if( params.fieldimage == "" ) {
                field_fits = Channel.fromPath("${params.publish_dir}/${params.label}/loadfits/field/${params.label}_field.fits")
            }
            else {
                field_fits = Channel.fromPath("${params.fieldimage}")
            }
        }
        else {
            empty_file = create_empty_file("file")
            field_fits = corr_field(
                "${label}_field", data, fcm, ra0, dec0, empty_file, polyco, 0, "N/A", "field"
            )
        }

        if ( params.calibrate ) {
            if( new File("${params.publish_dir}/${params.label}/finder/${params.label}.jmfit").exists()
                and new File("${params.publish_dir}/${params.label}/position/${params.label}_final_position.txt").exists() ) 
            {
                askap_frb_pos = Channel.fromPath("${params.publish_dir}/${params.label}/finder/${params.label}.jmfit")
            }
            else {
                if ( params.fieldflagfile == "" ) {
                    println "No field flag file!"
                    System.exit(1)
                }
                if ( params.flagfinder ) {
                    no_rfi_finder_fits = finder_fits
                }
                else {
                    no_rfi_finder_fits = subtract_rfi_finder(finder_fits, rfi_fits, subtractions, "finder")
                }

                bin_jmfits, bin_images, bin_regs, bin_mss = cal_finder(
                    no_rfi_finder_fits, flux_cal_solns, label, cpasspoly
                )

                askap_frb_pos = get_peak(
                    bin_jmfits.collect(), bin_images.collect(), 
                    bin_regs.collect(), bin_mss.collect()
                ).peak_jmfit

                field_sources = cal_field(
                    field_fits, flux_cal_solns, fieldflagfile, label, cpasspoly, askap_frb_pos
                ).jmfit

                apply_offset(field_sources, askap_frb_pos)
            }
        }

        if ( params.beamform ) {
            bform_frb(
                label, data, fcm, askap_frb_pos, flux_cal_solns, pol_cal_solns,
                num_ints, int_len, offset, dm, centre_freq, "-ds -t -XYIQUV"
            )
        }

    // emit:
        // true_pos = apply_offset.out
        // frb_htr_data = beamform_frb.out
}
