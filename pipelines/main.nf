nextflow.enable.dsl=2

include { process_flux_cal as fcal } from './process_flux_cal'
include { process_pol_cal as pcal } from './process_pol_cal'
include { process_frb as frb } from './process_frb'
include { create_empty_file } from './utils'

// Defaults
params.fluxflagfile = ""
params.polflagfile = ""
params.fieldflagfile = ""

params.calibrate = false
params.beamform = false

params.autoflag = true

params.nofrb = false        // can be convenient to not run frb processes
params.nopolcal = false     // some FRBs have no good pol cal

params.target = "FRB${params.label}"

params.out_dir = "${params.publish_dir}/${params.label}"

workflow {
    flux_cal_solns = fcal()

    if(params.nopolcal) {
        pol_cal_solns = create_empty_file("polcal.dat")
    }
    else {
        pol_cal_solns = pcal(
            flux_cal_solns,
        )
    }

    if(!params.nofrb) {
        frb(
            flux_cal_solns,
            pol_cal_solns,
        )
    }
}
