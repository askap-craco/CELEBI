//	This version is based on the branch developed by Paul Hancock
//	Recent updates are bering merged into it by AB

nextflow.enable.dsl=2
nextflow.enable.strict=true // be less generous

include { process_flux_cal as fcal1; process_flux_cal as fcal2 } from './process_flux_cal'
include { image_fluxcal } from './calibration'
include { process_pol_cal as pcal } from './process_pol_cal'
include { process_frb as frb } from './process_frb'
include { create_empty_file as empty1; create_empty_file as empty2; print_params as printpar;
    make_default_flags as makedefaultflgs; usable_ants as good_ants } from './utils'

utils_dir    = "${projectDir}/../utils/"
beamform_dir = "${projectDir}/../beamform/"
localise_dir = "${projectDir}/../localise/"

// Here's the main CELEBI workflow

workflow {

    if (params.askapbeam != "" && params.pols.size() != 1) {
        exit 1, "Error: --askapbeam is specified for single pol data, but params.pols does not have exactly 1 polarisation."
    }

    printpar()    

    makedefaultflgs()
    goodantfile = makedefaultflgs.out.goodants
    flgantfile = makedefaultflgs.out.flaggedants
    antcountfile = makedefaultflgs.out.antcounts
    antcnts = good_ants(goodantfile, antcountfile)
    antsfcal = antcnts.antsfcal.splitCsv(sep:",").flatten()
    antspcal = antcnts.antspcal.splitCsv(sep:",").flatten()
    antsfrb = antcnts.antsfrb.splitCsv(sep:",").flatten()
    antscomm = antcnts.antscomm.splitCsv(sep:",").flatten()

    if( params.localize || params.fcal ) {
        fcm_delayfix = fcal1(params.fcm).fcm_delayfix
        
        if(fcm_delayfix != "") {
            flux_cal_solns = fcal2(fcm_delayfix).flux_cal_solns
        }
    }
    else {
        fluxcal_solns_path = "${params.out_dir}/fluxcal/calibration_noxpol_${params.target}.tar.gz"
        fcm_delayfix_path = "${params.out_dir}/fluxcal/fcm_delayfix.txt"
        flux_cal_solns = file(fluxcal_solns_path)
        fcm_delayfix = file(fcm_delayfix_path)
    }

    if( params.localize || params.imgfcal ) {
        if ( params.noflag ) {
            fcalfits = file("${params.out_dir}/loadfits/fluxcal/${params.label}_fluxcal.fits")
        }
        else {
            fcalfits = file("${params.out_dir}/loadfits/fluxcal/${params.label}_fluxcal_f.fits")
        }
        
        image_fluxcal(
            fcalfits, flux_cal_solns, params.fluxflagfile
        )
    }
    
    if( params.nopolcal ) {
        pol_cal_solns = empty1("polcal.dat")
    }
    else if (params.psoln != "") {
        pol_cal_solns = file(params.psoln)
    }
    else {
        pol_cal_solns = pcal(
            flux_cal_solns, fcm_delayfix, antspcal
        )
    }

    if(!params.nofrb) {
        frb(
            flux_cal_solns,
            pol_cal_solns,
            fcm_delayfix,
            antsfrb
        )
    }
}
