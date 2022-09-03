nextflow.enable.dsl=2

include { process_flux_cal as fcal } from './process_flux_cal'
include { process_pol_cal as pcal } from './process_pol_cal'
include { process_frb as frb } from './process_frb'
include { create_empty_file } from './utils'

// Defaults
params.cpasspoly_fluxcal = 5
params.cpasspoly_polcal = 5
params.cpasspoly_frb = 5

params.fluxflagfile = ""
params.polflagfile = ""
params.fieldflagfile = ""

params.nocorrelate = false
params.calibrate = false
params.nocalibrate = false
params.beamform = false

params.nofrb = false        // can be convenient to not run frb processes
params.nopolcal = false     // some FRBs have no good pol cal

params.outdir = "${params.publish_dir}/${params.label}"

workflow {
    flux_cal_solns = fcal(
        "${params.label}_fluxcal",
        params.label,
        params.data_fluxcal,
        params.data_frb,
        params.snoopy,
        params.fcm,
        params.ra_fluxcal,
        params.dec_fluxcal,
        params.fluxflagfile,
        params.cpasspoly_fluxcal
    )
    if ( params.nopolcal ) {
        pol_cal_solns = create_empty_file("polcal.dat")
    }
    else {
        pol_cal_solns = pcal(
            "${params.label}_polcal",
            params.label,
            params.data_polcal,
            params.data_frb,
            params.snoopy,
            params.fcm,
            params.ra_polcal,
            params.dec_polcal,
            params.polflagfile,
            params.cpasspoly_polcal,
            flux_cal_solns,
            params.dm_polcal,
            params.centre_freq_polcal
        )
    }
    if(!params.nofrb){
        frb(
            params.label,
            params.data_frb,
            params.snoopy,
            params.fcm,
            params.ra_frb,
            params.dec_frb,
            params.fieldflagfile,
            flux_cal_solns,
            pol_cal_solns,
            params.cpasspoly_frb,
            params.dm_frb,
            params.centre_freq_frb
        )
    }
}
