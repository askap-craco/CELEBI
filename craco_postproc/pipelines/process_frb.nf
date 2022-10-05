nextflow.enable.dsl=2

include { create_empty_file } from './utils'
include { correlate as corr_finder; correlate as corr_rfi;
    correlate as corr_field; correlate as corr_htrgate; correlate as corr_htrrfi; 
    subtract_rfi as sub_rfi; subtract_rfi as sub_htrrfi } from './correlate'
include { image_finder; image_field; get_peak; image_htrgate } from './calibration'
include { find_offset; apply_offset; apply_offset as apply_offset_htr; 
    generate_binconfig } from './localise'
include { beamform as bform_frb; dedisperse; ifft; generate_dynspecs } from './beamform'

params.fieldimage = ""
params.flagfinder = ""
params.skiprfi = false

params.opt_DM = false
params.minDM = 0
params.maxDM = 10
params.DMstep = 0.01
params.opt_DM_dt = 100

params.opt_gate = false

beamform_dir = "$baseDir/../beamform/"

process plot {
    /*
        Plot dynamic spectra across different time resolutions to produce
        summary plot.

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            fnames_file: path
                File containing file names of dynamic spectra
            dynspecs: path
                Stokes parameter dynamic spectra
            centre_freq: val
                Central frequency of fine spectrum (MHz)
            dm: val
                Dispersion measure the data has been dedispersed to
            xy: path
                High-time resolution time series of X and Y pols
        
        Output:
            plot: path
                Plotted dynamic spectra across different time resolutions
            crops: path
                Directory containing cropped numpy files
    */
    publishDir "${params.publish_dir}/${params.label}/htr", mode: "copy"

    input:
        val label
        path fnames_file
        path dynspecs
        val centre_freq
        val dm
        path xy
    
    output:
        path "*.png"
        path "crops", emit: crops
        path "50us_crop_start_s.txt", emit: crop_start
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            module load gcc/9.2.0
            module load openmpi/4.0.2
            module load python/3.7.4
            module load numpy/1.18.2-python-3.7.4
            module load matplotlib/3.2.1-python-3.7.4
        fi
        args="-s $fnames_file"
        args="\$args -f $centre_freq"
        args="\$args -l $label"
        args="\$args -d $dm"
        args="\$args -x ${label}*X_t*npy"
        args="\$args -y ${label}*Y_t*npy"

        mkdir crops

        python3 $beamform_dir/plot.py \$args
        """
    
    stub:
        """
        touch stub.png
        mkdir crops
        touch 50us_crop_start_s.txt
        """
}

/*
process npy_to_archive {
    /*
        Convert X and Y numpy time series to a PSRCHIVE archive

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            crops: path
                Directory containing cropped numpy files to convert
            startmjd: val
                Earliest data start time in Modified Julian Day (MJD) 
            centre_freq: val
                Central frequency of fine spectrum (MHz)
            final_position: path
                Text file containing final position and error
            
        Output:
            fils: path
                Directory containing converted filterbank files
    
    publishDir "${params.publish_dir}/${params.label}/htr/fils", mode: "copy"

    input:
        val label
        path crops
        val startmjd
        val centre_freq
        path final_position
    
    output:
        path "*fil"
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            module load anaconda3/5.1.0
            source activate $launchDir/envs/psrchive
        fi

        #convert_addpol
        #fill_header
        #cat outputs of ^ into outfile

        #dspsr
        #pam
        """

    stub:
        """
        touch stub.fil
        """
}
*/

process find_DM_opt {
    /*
        Optimise DM for S/N. Works under the assumption that the current DM
        is an underestimate

        Input
            crops: path
                Cropped FRB data
            dm: val
                Current DM

        Output
            stdout
                S/N maximising DM
            opt_DM_plot
                max(I) vs DM plot
    */
    publishDir "${params.publish_dir}/${params.label}/htr", mode: "copy"

    input:
        path crops
        val dm

    output:
        env dmopt, emit: dm_opt
        path "*png"

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            module load gcc/9.2.0
            module load openmpi/4.0.2
            module load python/3.7.4
            module load numpy/1.18.2-python-3.7.4
            module load matplotlib/3.2.1-python-3.7.4
        fi

        args="-x $crops/${params.label}_${dm}_X.npy"
        args="\$args -y $crops/${params.label}_${dm}_Y.npy"
        args="\$args -d $params.minDM"
        args="\$args -D $params.maxDM"
        args="\$args -s $params.DMstep"
        args="\$args --DM0 $dm"
        args="\$args --f0 $params.centre_freq_frb"
        args="\$args --dt $params.opt_DM_dt"

        dmopt=`python3 $beamform_dir/opt_DM.py \$args`
        """
    
    stub:
        """
        dmopt=$dm
        touch stub.png
        """
}

workflow optimise_DM {
    /*
        After initial beamforming, optimise the DM and re-generate plots
    */
    take:
        pre_dedisp
        crops
        pol_cal_solns
        ds_args
    
    main:
        find_DM_opt(crops, params.dm_frb)
        dm_opt = find_DM_opt.out.dm_opt
        dedisperse(
            params.label, dm_opt, params.centre_freq_frb, pre_dedisp
        )
        ifft(params.label, dedisperse.out, dm_opt)
        xy = ifft.out.collect()
        generate_dynspecs(
            params.label, xy, ds_args, params.centre_freq_frb, dm_opt, pol_cal_solns
        )
        plot(
            params.label, generate_dynspecs.out.dynspec_fnames, 
            generate_dynspecs.out.data, params.centre_freq_frb, dm_opt,
            xy
        )
    
    emit:
        dm_opt
        crops = plot.out.crops
        crop_start = plot.out.crop_start
}

process mjd_prof {
    /*
        Create profile as function of MJD

        Input
            crops: path
                Crops directory from plot
            crop_start: path
                File containing start time of 50us crop relative to full data 
                in seconds
            
        Output
            prof: path
                Two-column space separated file containing MJD and 50us profile
                respectively
    */
    input:
        path crops
        path crop_start
    
    output:
        path "prof.txt", emit: prof

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            module load gcc/9.2.0
            module load openmpi/4.0.2
            module load python/3.7.4
            module load numpy/1.18.2-python-3.7.4
        fi

        python3 $beamform_dir/mjd_prof.py $params.data_frb $crops/*_50us_I.npy $crop_start
        """
    script:
        """
        touch prof.txt
        """
}

process update_polyco {
    /*
        Edit a polyco file to replace the DM with a new value

        Input
            polyco: path
                Polyco file to edit
            dm: val
                DM value to insert into polyco
        
        Output
            craftfrb.polyco: path
                Edited polyco file
    */
    input:
        path polyco, stageAs: "old.polyco"
        val dm
    
    output:
        path "craftfrb.polyco"
    
    script:
        """
        head -1 $polyco | awk '\$5="$dm"' > craftfrb.polyco
        tail -2 $polyco >> craftfrb.polyco
        """
    
    stub:
        """
        cp old.polyco craftfrb.polyco
        """
}

process htr_to_binconfig {
    /*
        Write a binconfig file with a matched filter for a provided profile
        as a function of MJD

        Input
            prof: path
                High time resolution profile as function of MJD
            polyco: path
                Polyco file
        
        Output:
            htr gate binconfig: path
                Binconfig containing matched filter for high time res gate
    */
    input:
        path prof
        path polyco
    
    output:
        path "craftfrb.htrgate.binconfig", emit: htrgate
        path "craftfrb.htrrfi.binconfig", emit: htrrfi
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            module load gcc/9.2.0
            module load openmpi/4.0.2
            module load python/3.7.4
            module load numpy/1.18.2-python-3.7.4
            module load matplotlib/3.2.1-python-3.7.4
        fi

        python3 $beamform_dir/htr2binconfig.py $prof $polyco
        """
    
    stub:
        """
        touch craftfrb.htrgate.binconfig
        touch craftfrb.htrrfi.binconfig
        """
}

workflow optimise_gate {
    /*
        Create an optimised matched filter binconfig based on beamformed high
        time resolution data and optimised DM
    */
    take:
        crops
        crop_start
        polyco
        dm
    
    main:
        new_polyco = update_polyco(polyco, dm)
        prof = mjd_prof(crops, crop_start)
        htr_to_binconfig(prof, polyco)
    
    emit:
        htrgate = htr_to_binconfig.out.htrgate
        htrrfi = htr_to_binconfig.out.htrrfi
        polyco = new_polyco
}

workflow process_frb {
    /*
        Process voltages to obtain an FRB position

        Take
            flux_cal_solns: path
                Flux calibrator solutions tarball
            pol_cal_solns: path
                Polarisation calibration solutions in a text file
    */
    take:
        flux_cal_solns
        pol_cal_solns

    main:
        binconfig = generate_binconfig()
        empty_file = create_empty_file("file")

        if(!params.opt_gate){    
            // Correlate finder
            finder_fits_path = "${params.publish_dir}/${params.label}/loadfits/finder/finderbin20.fits"
            if(new File(finder_fits_path).exists()) {
                finder_fits = Channel.fromPath(
                    "${params.publish_dir}/${params.label}/loadfits/finder/finderbin*.fits"
                )
            }
            else {
                finder_fits = corr_finder(
                    "finder", params.data_frb, params.ra_frb, params.dec_frb, 
                    binconfig.finder, binconfig.polyco, binconfig.int_time, "finder"
                )
            }

            // Correlate RFI (if not directly flagging finder)
            rfi_fits_path = "${params.publish_dir}/${params.label}/loadfits/rfi/${params.label}_rfi.fits"
            if ( new File(rfi_fits_path).exists() ) {
                rfi_fits = Channel.fromPath(rfi_fits_path)
            }
            else {
                if(!params.skiprfi) {
                    rfi_fits = corr_rfi(
                        "${params.label}_rfi", params.data_frb, params.ra_frb, 
                        params.dec_frb, binconfig.rfi, binconfig.polyco, binconfig.int_time, "rfi"
                    )
                }
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
            field_fits = corr_field(
                "${params.label}_field", params.data_frb, params.ra_frb, 
                params.dec_frb, empty_file, empty_file, 0, "field"
            )
        }

        // Calibrate (i.e. image finder and field)
        frb_jmfit_path = "${params.publish_dir}/${params.label}/finder/${params.label}.jmfit"
        offset_path = "${params.publish_dir}/${params.label}/position/offset0.dat"
        frb_pos_path = "${params.publish_dir}/${params.label}/position/${params.label}_final_position.txt"
        if(new File(frb_jmfit_path).exists()) {
            askap_frb_pos = Channel.fromPath(frb_jmfit_path)
        }
        else if(new File(frb_pos_path).exists()) {
            final_position = Channel.fromPath(frb_pos_path)
        }
        if(params.calibrate) {
            if((params.fieldimage == "") && (params.fieldflagfile == "")){
                println "No field flag file!"
                System.exit(1)
            }
            
            if(!params.opt_gate){
                if(params.skiprfi){
                    no_rfi_finder_fits = finder_fits
                }
                else {
                    no_rfi_finder_fits = sub_rfi(
                        finder_fits, rfi_fits, binconfig.subtractions
                    )                
                }

                bins_out = image_finder(
                    no_rfi_finder_fits, flux_cal_solns
                )
                bin_jmfits = bins_out.jmfit
                bin_fits_images = bins_out.fits_image
                bin_regs = bins_out.reg
                bin_mss = bins_out.ms

                askap_frb_pos = get_peak(
                    bin_jmfits.collect(), bin_fits_images.collect(), 
                    bin_regs.collect(), bin_mss.collect()
                ).peak_jmfit
            }
            else {
                askap_frb_pos = empty_file
            }

            field_sources = image_field(
                field_fits, flux_cal_solns, params.fieldflagfile, askap_frb_pos
            ).jmfit

            offset = find_offset(field_sources).offset

            if(!params.opt_gate){
                final_position = apply_offset(offset, askap_frb_pos).final_position
            }
        }
        else if(new File(offset_path).exists()) {
            offset = Channel.fromPath(offset_path)
        }

        if(params.beamform) {
            bform_frb(
                params.label, params.data_frb, askap_frb_pos, flux_cal_solns, 
                pol_cal_solns, params.dm_frb, params.centre_freq_frb, "-ds -t -XYIQUV"
            )
            plot(
                params.label, bform_frb.out.dynspec_fnames, bform_frb.out.htr_data,
                params.centre_freq_frb, params.dm_frb, bform_frb.out.xy
            )

            if(params.opt_DM) {
                optimise_DM(
                    bform_frb.out.pre_dedisp, plot.out.crops, pol_cal_solns, 
                    "-ds -t -XYIQUV"
                )
                dm = optimise_DM.out.dm_opt
                crops = optimise_DM.out.crops
                crop_start = optimise_DM.out.crop_start
            }
            else {
                dm = params.dm_frb
                crops = plot.out.crops
                crop_start = plot.out.crop_start
            }

            if(params.opt_gate) {
                opt_gate = optimise_gate(crops, crop_start, binconfig.polyco, dm)
                htrgate_fits = corr_htrgate(
                    "${params.label}_htrgate", params.data_frb, params.ra_frb, 
                    params.dec_frb, opt_gate.htrgate, opt_gate.polyco, 
                    binconfig.int_time, "htrgate"
                )
                if(!params.skiprfi) {
                    htrrfi_fits = corr_htrrfi(
                        "${params.label}_htrrfi", params.data_frb, params.ra_frb, 
                        params.dec_frb, opt_gate.htrrfi, opt_gate.polyco, 
                        binconfig.int_time, "htrrfi"
                    )
                }
                if(params.skiprfi){
                    no_rfi_htrgate_fits = htrgate_fits
                }
                else {
                    no_rfi_htrgate_fits = sub_htrrfi(    
                        htrgate_fits, htrrfi_fits, empty_file
                    )                
                }

                image_htrgate(
                    no_rfi_htrgate_fits, flux_cal_solns
                )

                final_position = apply_offset_htr(
                    offset, image_htrgate.out.jmfit
                )
            }

            //npy2fil(params.label, plot.out.crops, 0, params.centre_freq_frb, final_position)
        }
}
