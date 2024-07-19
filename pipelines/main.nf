nextflow.enable.dsl=2

include { process_flux_cal as fcal1; process_flux_cal as fcal2 } from './process_flux_cal'
include { process_pol_cal as pcal } from './process_pol_cal'
include { process_frb as frb } from './process_frb'
include { create_empty_file as empty1; create_empty_file as empty2 } from './utils'

// Defaults
params.fluxflagfile = ""
params.polflagfile = ""
params.fieldflagfile = ""
params.calibrate = false
params.beamform = false
params.noflag = false       // don't automatically flag
params.nofrb = false        // can be convenient to not run frb processes
params.nopolcal = false     // some FRBs have no good pol cal
params.target = "FRB${params.label}"
params.out_dir = "${params.publish_dir}/${params.label}"
params.psoln = ""
params.nbits = 4
params.numfinderbins = 7
params.searchms = 70
params.uselocalcatalog = false
params.referencecatalog = "RACS"
params.nants = 2
params.nants_fcal = params.nants

workflow {
    fcm_delayfix = fcal1(params.fcm).fcm_delayfix
    //fcm_delayfix = fcal1("/fred/oz313/data/frb210912/fcm.txt.32063").fcm_delayfix
    //println "TESTING"
    println fcm_delayfix
    println params.fcm
    if(fcm_delayfix != "") {
        flux_cal_solns = fcal2(fcm_delayfix).flux_cal_solns
    }

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
