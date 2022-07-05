nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_pcal } from './correlate'
include { beamform as bform_pcal } from './beamform'
include { apply_flux_cal_solns_polcal as apply_cal_pcal; determine_pol_cal_solns as get_cal_pcal } from './calibration'
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
        empty_file = create_empty_file("file")
        if ( new File("${params.publish_dir}/${params.label}/loadfits/polcal/${params.label}_polcal.fits").exists() ) {
            fits = Channel.fromPath("${params.publish_dir}/${params.label}/loadfits/polcal/${params.label}_polcal.fits")
        }
        else {
            binconfig = generate_binconfig(data_frb, snoopy)
            fits = corr_pcal(
                label, data, fcm, ra0, dec0, empty_file, binconfig.polyco, 0, polflagfile, "polcal"
            )
        }

        if ( params.calibrate ) {
            if ( new File("${params.publish_dir}/${params.label}/polcal/polcal.jmfit").exists() ) {
                pos = Channel.fromPath("${params.publish_dir}/${params.label}/polcal/polcal.jmfit")
            }
            else {
                if ( params.polflagfile == "" ) {
                    println "No polcal flag file!"
                    System.exit(1)
                }
                pos = apply_cal_pcal(
                    fits, flux_cal_solns, polflagfile, target, cpasspoly
                ).jmfit
            }
        }

        // if ( params.beamform ) {
        //     htr_data = bform_pcal(
        //         label, data, fcm, pos, flux_cal_solns, empty_file,
        //         num_ints, int_len, offset, dm, centre_freq, "-ds -IQUV"
        //     )
        //     pol_cal_solns = get_cal_pcal(htr_data).pol_cal_solns
        // }   
        // else {
        pol_cal_solns = ""
        // }

    emit:
        pol_cal_solns
}
