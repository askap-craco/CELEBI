nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_fcal } from './correlate'
include { determine_flux_cal_solns as cal_fcal } from './calibration'
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
        if( new File("${params.publish_dir}/${params.label}/loadfits/fluxcal/${params.label}_fluxcal.fits").exists() ) {
            fits = Channel.fromPath("${params.publish_dir}/${params.label}/loadfits/fluxcal/${params.label}_fluxcal.fits")
        }
        else {
            binconfig = generate_binconfig(data_frb, snoopy)
            empty_binconfig = create_empty_file("binconfig")
            fits = corr_fcal(
                label, data, fcm, ra, dec, empty_binconfig, binconfig.polyco, 0, 
                fluxflagfile, "fluxcal"
            )
        }

        if ( params.calibrate ) {
            if( new File("${params.publish_dir}/${params.label}/fluxcal/calibration_noxpol_${target}.tar.gz").exists() ) {
                flux_cal_solns = Channel.fromPath("${params.publish_dir}/${params.label}/fluxcal/calibration_noxpol_${target}.tar.gz")
            }
            else {
                if ( params.fluxflagfile == "" ) {
                    println "No fluxcal flag file!"
                    System.exit(1)
                }
                flux_cal_solns = cal_fcal(fits, fluxflagfile, target, cpasspoly).solns
            }
        }
        else {
            flux_cal_solns = ""
        }
    emit:
        flux_cal_solns
}
