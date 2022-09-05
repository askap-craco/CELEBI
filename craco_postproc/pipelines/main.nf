nextflow.enable.dsl=2

include { process_flux_cal as fcal } from './process_flux_cal'
include { process_pol_cal as pcal } from './process_pol_cal'
include { process_frb as frb } from './process_frb'
include { create_empty_file } from './utils'
include { generate_binconfig } from './localise'

// Defaults
params.fluxflagfile = ""
params.polflagfile = ""
params.fieldflagfile = ""

params.nocorrelate = false
params.calibrate = false
params.nocalibrate = false
params.beamform = false

params.nofrb = false        // can be convenient to not run frb processes
params.nopolcal = false     // some FRBs have no good pol cal

params.target = "FRB${params.label}"

params.outdir = "${params.publish_dir}/${params.label}"

workflow {
    binconfig = generate_binconfig()

    flux_cal_solns = fcal(
        binconfig.polyco,
    )

    if(params.nopolcal) {
        pol_cal_solns = create_empty_file("polcal.dat")
    }
    else {
        pol_cal_solns = pcal(
            binconfig.polyco,
            flux_cal_solns,
        )
    }

    if(!params.nofrb) {
        frb(
            binconfig,              // expands into 6 Channels
            flux_cal_solns,
            pol_cal_solns,
        )
    }
}
