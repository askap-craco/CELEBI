nextflow.enable.dsl=2

include { process_flux_cal } from './process_flux_cal'
include { process_pol_cal } from './process_pol_cal'
include { process_frb } from './process_frb'
include { create_empty_file } from './utils'

// Defaults
params.cpasspoly_fluxcal = 5
params.cpasspoly_polcal = 5
params.cpasspoly_frb = 5

params.num_ints_polcal = 1
params.int_len_polcal = 44000
params.offset_polcal = 0

params.num_ints_frb = 1
params.int_len_frb = 44000
params.offset_frb = 0

params.fluxflagfile = ""
params.polflagfile = ""
params.fieldflagfile = ""

params.nocorrelate = false
params.calibrate = false
params.beamform = false

workflow {
    flux_cal_solns = process_flux_cal(
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
}
