nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_finder; correlate as corr_rfi;
    correlate as corr_field; subtract_rfi } from './correlate'
include { apply_flux_cal_solns_finder as cal_finder;
    apply_flux_cal_solns_field as cal_field; get_peak } from './calibration'
include { apply_offset; generate_binconfig } from './localise'
include { beamform as bform_frb } from './beamform'

params.fieldimage = ""
params.flagfinder = false

workflow process_frb {
    /*
        Process voltages to obtain an FRB position

        Take
            label: val
                FRB name and context as a string (no spaces)
            data: val
                Absolute path to FRB data base directory (the dir. with the ak*
                directories)
            snoopy: val
                Absolute path to snoopyv2.log file of FRB trigger
            fcm: val
                Absolute path to fcm (hardware delays) file
            ra0: val
                FRB right ascension initial guess as "hh:mm:ss"
            dec0: val
                FRB declination initial guess as "dd:mm:ss"
            fieldflagfile: val
                Absolute path to AIPS flag file for field visibilities. If set 
                to a blank string, the workflow will end before calibrating.
            flux_cal_solns: path
                Flux calibrator solutions tarball
            pol_cal_solns: path
                Polarisation calibration solutions in a text file
            cpasspoly: val
                Order of polynomial to fit bandpass with
            num_ints: val
                TODO: deprecate
            int_len: val
                TODO: deprecate
            offset: val
                TODO: deprecate
            dm: val
                DM to dedisperse to in pc/cm3
            centre_freq: val
                Central frequency of data in MHz
    */
    take:
        label
        data
        snoopy
        fcm
        ra0
        dec0
        fieldflagfile
        flux_cal_solns
        pol_cal_solns
        cpasspoly
        num_ints
        int_len
        offset
        dm
        centre_freq

    main:
        binconfig = generate_binconfig(data, snoopy)
        binconfig_finder = binconfig.finder
        binconfig_rfi = binconfig.rfi
        subtractions = binconfig.subtractions
        polyco = binconfig.polyco
        int_time = binconfig.int_time

        // Correlate finder
        finder_fits_path = "${params.publish_dir}/${params.label}/loadfits/finder/finderbin20.fits"
        if(new File(finder_fits_path).exists()) {
            finder_fits = Channel.fromPath(finder_fits_path)
        }
        else {
            finder_fits = corr_finder(
                "finder", data, fcm, ra0, dec0, binconfig_finder, polyco, int_time, 
                "finder"
            )
        }

        // Correlate RFI (if not directly flagging finder)
        rfi_fits_path = "${params.publish_dir}/${params.label}/loadfits/rfi/${params.label}_rfi.fits"
        if ( new File(rfi_fits_path).exists() ) {
            rfi_fits = Channel.fromPath(rfi_fits_path)
        }
        else {
            if(!params.flagfinder) {
                rfi_fits = corr_rfi(
                    "${label}_rfi", data, fcm, ra0, dec0, binconfig_rfi, polyco, 
                    int_time, "rfi"
                )
            }
        }

        // Correlate field (if not using deep field image)
        field_fits_path = "${params.publish_dir}/${params.label}/loadfits/field/${params.label}_field.fits"
        if((params.fieldimage != "") or new File(field_fits_path).exists() ) {
            if(params.fieldimage == "") {
                field_fits = Channel.fromPath(field_fits_path)
            }
            else {
                field_fits = Channel.fromPath("${params.fieldimage}")
            }
        }
        else {
            empty_file = create_empty_file("file")
            field_fits = corr_field(
                "${label}_field", data, fcm, ra0, dec0, empty_file, polyco, 0, "field"
            )
        }

        // Calibrate (i.e. image finder and field)
        if(params.calibrate) {
            frb_jmfit_path = "${params.publish_dir}/${params.label}/finder/${params.label}.jmfit"
            frb_pos_path = "${params.publish_dir}/${params.label}/position/${params.label}_final_position.txt"
            if(new File(frb_jmfit_path).exists() && new File(frb_pos_path).exists()) {
                askap_frb_pos = Channel.fromPath(frb_jmfit_path)
            }
            else {
                if((params.fieldimage == "") && (params.fieldflagfile == "")){
                    println "No field flag file!"
                    System.exit(1)
                }
                
                if(params.flagfinder){
                    no_rfi_finder_fits = finder_fits
                }
                else {
                    no_rfi_finder_fits = subtract_rfi(
                        finder_fits, rfi_fits, subtractions
                    )                
                }

                bins_out = cal_finder(
                    no_rfi_finder_fits, flux_cal_solns, label, cpasspoly
                )
                bin_jmfits = bins_out.jmfit
                bin_fits_images = bins_out.fits_image
                bin_regs = bins_out.reg
                bin_mss = bins_out.ms

                askap_frb_pos = get_peak(
                    bin_jmfits.collect(), bin_fits_images.collect(), 
                    bin_regs.collect(), bin_mss.collect()
                ).peak_jmfit

                field_sources = cal_field(
                    field_fits, flux_cal_solns, fieldflagfile, label, cpasspoly, 
                    askap_frb_pos
                ).jmfit

                apply_offset(field_sources, askap_frb_pos)
            }
        }

        if(params.beamform) {
            bform_frb(
                label, data, fcm, askap_frb_pos, flux_cal_solns, pol_cal_solns,
                num_ints, int_len, offset, dm, centre_freq, "-ds -t -XYIQUV"
            )
        }
}
