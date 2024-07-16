nextflow.enable.dsl=2
nextflow.enable.strict=true // be less generous

include { process_flux_cal as fcal1; process_flux_cal as fcal2 } from './process_flux_cal'
include { process_pol_cal as pcal } from './process_pol_cal'
include { process_frb as frb } from './process_frb'
include { create_empty_file as empty1; create_empty_file as empty2 } from './utils'

// Defaults
params.fluxflagfile = ""
params.polflagfile = ""
params.fieldflagfile = ""
params.calibrate = true
params.beamform = true
params.noflag = false       // don't automatically flag
params.nofrb = false        // can be convenient to not run frb processes
params.nopolcal = false     // some FRBs have no good pol cal
params.target = "FRB${params.label}"
params.out_dir = "${params.publish_dir}/${params.label}"
params.psoln = ""
params.nants = 2
params.nants_fcal = params.nants

workflow {
    fcm_delayfix = fcal1(params.fcm).fcm_delayfix

    // the following will always fail if calibrate=false, why is thise even an option?
    if(fcm_delayfix != "") {
        flux_cal_solns = fcal2(fcm_delayfix).flux_cal_solns
    }
    // else there is no fcm_delayfix and frb will fail

    if(params.nopolcal) {
        pol_cal_solns = empty1("polcal.dat")
    }
    else if (params.psoln != "") {
        pol_cal_solns = Channel.fromPath(params.psoln)
    }
    else {
        pol_cal_solns = pcal(
            flux_cal_solns, fcm_delayfix
        )
    }

    if(!params.nofrb) {
        frb(
            flux_cal_solns,
            pol_cal_solns,
            fcm_delayfix
        )
    }
}
