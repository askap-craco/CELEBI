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

params.image_all_bins = false

params.ICS_DMrange = 100
params.ICS_DMstep = 0.1

params.opt_DM = false
params.minDM = 0
params.maxDM = 10
params.DMstep = 0.01
params.opt_DM_dt = 100

params.opt_gate = false

params.pols = ['X', 'Y']
polarisations = Channel
    .fromList(params.pols)
params.nants = 2
antennas = Channel
    .of(0..params.nants-1)

//beamform_dir = "$baseDir/../beamform/"
beamform_dir = "/fred/oz002/askap/craft/craco/CELEBI/craco_postproc/beamform"
localise_dir = "/fred/oz002/askap/craft/craco/CELEBI/craco_postproc/localise"

process load_coarse_dynspec {
    /*
        Incoherently create a 1 ms dynamic spectrum from voltages for a given
        polarisation and antenna

        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            data: val
                Absolute path to data base directory (the dir. with the ak* 
                directories)
            pol: val
                One of "X" or "Y" for the current polarisation being beamformed
            ant_idx: val
                Zero-based index of the antenna being beamformed

        Output
            data: path
                1 ms time resolution dynamic spectrum
            time: path
                Time axis in MJD
    */

    input:
        val label
        val data
        each pol
        each ant_idx

    output:
        path "${label}_ICS_${pol}*${ant_idx}.npy", emit: data
        path "t_mjd.npy", emit: time

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi

        # create calcfiles for antenna delay
        startmjd=`python3 /fred/oz002/askap/craft/craco/CELEBI/craco_postproc/localise/get_start_mjd.py $data`

        export CRAFTCATDIR="."

        # Run processTimeStep.py with the --calconly flag to stop once calcfile is 
        # written
        args="-t $data"
        args="\$args --ra $params.ra_frb"
        args="\$args -d$params.dec_frb"
        args="\$args -f $params.fcm"
        args="\$args -b 4"
        args="\$args --card 1"
        args="\$args -k"
        args="\$args --name=210117_ICS"
        args="\$args -o ."
        args="\$args --freqlabel c1_f0"
        args="\$args --dir=/fred/oz002/askap/craft/craco/CELEBI/craco_postproc/difx"
        args="\$args --calconly"
        args="\$args --startmjd \$startmjd"

        echo "python3 /fred/oz002/askap/craft/craco/CELEBI/craco_postproc/localise/processTimeStep.py \$args"
        python3 /fred/oz002/askap/craft/craco/CELEBI/craco_postproc/localise/processTimeStep.py \$args


        mkdir delays
        args="-d $data"
        args="\$args --parset $params.fcm"
        args="\$args --calcfile c1_f0/craftfrb.im"
        args="\$args -o ${label}_ICS"
        args="\$args --ics"
        args="\$args --cpus=8"
        args="\$args --pol=$pol"
        args="\$args --an=$ant_idx"

        echo "python3 $beamform_dir/craftcor_tab.py \$args"
        python3 $beamform_dir/craftcor_tab.py \$args
        """
    
    stub:
        """
        touch ${label}_ICS_${pol}_${ant_idx}.npy
        touch t_mjd.npy
        """
}

process refine_candidate {
    /*
        Sum incoherent dynamic spectra, search for FRB, and refine snoopy
        candidate
    
        Input
            label: val
                FRB name and context of process instance as a string (no
                spaces)
            ics_dynspecs: path
                All the incoherent dynamic spectra
            t_mjd: path
                ICS dynamic spectra time axis in MJD
            snoopy: path
                Initial detection snoopy candidate
    */
    publishDir "${params.publish_dir}/${params.label}/ics"

    input:
        val label
        path ics_dynspecs
        path t_mjd
        path snoopy

    output:
        path "${label}_ICS.npy", emit: sum_ics
        path "${label}.cand", emit: cand
        path "*png", emit: plots

    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi

        python3 $localise_dir/sum_ics.py ${label}_ICS.npy ${label}_ICS_*.npy

        args="--ds ${label}_ICS.npy"
        args="\$args -s $snoopy"
        args="\$args -t $t_mjd"
        args="\$args -f $params.centre_freq_frb"
        args="\$args --DMrange=$params.ICS_DMrange"
        args="\$args --DMstep=$params.ICS_DMstep"
        args="\$args -o ${label}.cand"

        python3 $localise_dir/search_ics.py \$args
        """

    stub:
        """
        touch ${label}_ICS.npy
        touch ${label}.cand
        touch stub.png
        """
}

process get_beam_centre {
    /*
        Parse VCRAFT headers to get the beam centre

        Output
            ra: env
                Beam centre right ascension (hms)
            dec: env
                Beam centre declination (dms)
    */
    output:
        env ra, emit: ra
        env dec, emit: dec
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi

        # find a header file
        ant_pattern="${params.data_frb}/ak*"
        ants=( \$ant_pattern )
        first_ant=`echo \$ants`
        beam_pattern="\$first_ant/beam*"
        beams=( \$beam_pattern )
        first_beam=`echo \$beams`
        header=`ls \$first_beam/*c1_f0*hdr`

        # beam centre in degrees
        ra_beam_deg=`grep BEAM_RA \$header | cut -d " " -f 2`
        dec_beam_deg=`grep BEAM_DEC \$header | cut -d " " -f 2`

        radec_beam=`python $localise_dir/get_beam_radec.py \$ra_beam_deg \$dec_beam_deg`
        ra=`echo \$radec_beam | cut -d " " -f 1 | tr h : | tr m : | tr s 0`
        dec=`echo \$radec_beam | cut -d " " -f 2 | tr d : | tr m : | tr s 0`
        """
    
    stub:
        """
        ra="00:00:00"
        dec="00:00:00"
        """
}

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
            time: path
                Time in MJD (1 ms steps) from coarse dynspecs
            cand: path
                Refined candidate for FRB from ICS search
        
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
        path time
        path cand
    
    output:
        path "*.png"
        path "crops", emit: crops
        path "50us_crop_start_s.txt", emit: crop_start
        path "crops/*50us_I.npy", emit: crop_50us
    
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
        args="\$args -t $time"
        args="\$args -c $cand"

        mkdir crops

        python3 $beamform_dir/plot.py \$args
        """
    
    stub:
        """
        touch stub.png
        mkdir crops
        touch crops/stub_50us_I.npy
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
        crop_50us = plot.out.crop_50us
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
        path crop_50us
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

        python3 $beamform_dir/mjd_prof.py $params.data_frb $crop_50us $crop_start
        """

    stub:
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
        head -2 $polyco | tail -1 | awk '\$6="1104.000"' >> craftfrb.polyco
        tail -1 $polyco >> craftfrb.polyco
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
        crop_50us
        crop_start
        polyco
        dm
    
    main:
        new_polyco = update_polyco(polyco, dm)
        prof = mjd_prof(crop_50us, crop_start)
        htr_to_binconfig(prof, new_polyco)
    
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
        refined_candidate_path = "${params.publish_dir}/${params.label}/ics/${params.label}.cand"
        if ( new File(refined_candidate_path).exists()) {
            refined_candidate = Channel.fromPath(refined_candidate_path)
        }
        else {
            coarse_ds = load_coarse_dynspec(params.label, params.data_frb, polarisations, 
                                            antennas)
            refine_candidate(params.label, coarse_ds.data.collect(), 
                             coarse_ds.time.first(), params.snoopy)
            refined_candidate = refine_candidate.out.cand
        }
        binconfig = generate_binconfig(refined_candidate)
        empty_file = create_empty_file("file")

        if(!params.opt_gate){    
            // Correlate finder
            finder_fits_path = "${params.publish_dir}/${params.label}/loadfits/finder/finderbin07.fits"
            if(new File(finder_fits_path).exists()) {
                finder_fits = Channel.fromPath(
                    "${params.publish_dir}/${params.label}/loadfits/finder/finderbin*.fits"
                )
                centre_bin_fits = Channel.fromPath(
                    "${params.publish_dir}/${params.label}/loadfits/finder/finderbin04.fits"
                )
            }
            else {
                (finder_fits, centre_bin_fits) = corr_finder(
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
                    ).fits
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
            beam_centre = get_beam_centre()
            field_fits = corr_field(
                "${params.label}_field", params.data_frb, beam_centre.ra, 
                beam_centre.dec, empty_file, empty_file, empty_file, "field"
            ).fits
        }

        // Calibrate (i.e. image finder and field)
        frb_jmfit_path = "${params.publish_dir}/${params.label}/finder/${params.label}.jmfit"
        offset_path = "${params.publish_dir}/${params.label}/position/offset0.dat"
        frb_pos_path = "${params.publish_dir}/${params.label}/position/${params.label}_final_position.txt"
        if(new File(frb_jmfit_path).exists()) {
            askap_frb_pos = Channel.fromPath(frb_jmfit_path)
        }
        if(new File(frb_pos_path).exists()) {
            final_position = Channel.fromPath(frb_pos_path)
        }
        if(params.calibrate) {
            if((params.fieldimage == "") && (params.fieldflagfile == "")){
                println "No field flag file!"
                System.exit(1)
            }
            
            if(!params.opt_gate){
                if(params.image_all_bins) {
                    bins_to_image = finder_fits
                }
                else {
                    bins_to_image = centre_bin_fits
                }

                if(params.skiprfi){
                    no_rfi_finder_fits = bins_to_image
                }
                else {
                    no_rfi_finder_fits = sub_rfi(
                        bins_to_image, rfi_fits, binconfig.subtractions
                    )                
                }

                bins_out = image_finder(
                    no_rfi_finder_fits, flux_cal_solns
                )
                bin_jmfits = bins_out.jmfit
                bin_fits_images = bins_out.fits_image
                bin_regs = bins_out.reg
                bin_mss = bins_out.ms

                if(params.image_all_bins) {
                    askap_frb_pos = get_peak(
                        bin_jmfits.collect(), bin_fits_images.collect(), 
                        bin_regs.collect(), bin_mss.collect()
                    ).peak_jmfit
                }
                else {
                    askap_frb_pos = bin_jmfits
                }
            }
            if(new File(frb_jmfit_path).exists()) {
                askap_frb_pos = Channel.fromPath(frb_jmfit_path)
            }
            // else {
            //     askap_frb_pos = empty_file
            // }

            if(new File(offset_path).exists()) {
                offset = Channel.fromPath(offset_path)
            }
            else {
                field_sources = image_field(
                    field_fits, flux_cal_solns, params.fieldflagfile, askap_frb_pos
                ).jmfit

                offset = find_offset(field_sources).offset
            }

            if(!params.opt_gate){
                final_position = apply_offset(offset, askap_frb_pos).final_position
            }
        }
        else if(new File(offset_path).exists()) {
            offset = Channel.fromPath(offset_path)
        }

        if(params.beamform) {
            crop_50us_path = "${params.publish_dir}/${params.label}/htr/crops/${params.label}_${params.dm_frb}_50us_I.npy"
            crop_start_path = "${params.publish_dir}/${params.label}/htr/50us_crop_start_s.txt"
            if(new File(crop_50us_path).exists() && new File(crop_start_path).exists()) {
                crop_50us = Channel.fromPath(crop_50us_path)
                crop_start = Channel.fromPath(crop_start_path)
            }
            else {
                bform_frb(
                    params.label, params.data_frb, askap_frb_pos, flux_cal_solns, 
                    pol_cal_solns, params.dm_frb, params.centre_freq_frb, "-ds -t -XYIQUV"
                )
                plot(
                    params.label, bform_frb.out.dynspec_fnames, bform_frb.out.htr_data,
                    params.centre_freq_frb, params.dm_frb, bform_frb.out.xy,
                    coarse_ds.time.first(), refined_candidate
                )
                crops = plot.out.crops
                crop_start = plot.out.crop_start
                crop_50us = plot.out.crop_50us
            }

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
            }

            if(params.opt_gate) {
                opt_gate = optimise_gate(crop_50us, crop_start, binconfig.polyco, dm)
                htrgate_fits_path = "${params.publish_dir}/${params.label}/loadfits/htrgate/${params.label}_htrgate.fits"
                if(new File(htrgate_fits_path).exists()) {
                    htrgate_fits = Channel.fromPath(htrgate_fits_path)
                }
                else {
                    htrgate_fits = corr_htrgate(
                        "${params.label}_htrgate", params.data_frb, params.ra_frb, 
                        params.dec_frb, opt_gate.htrgate, opt_gate.polyco, 
                        binconfig.int_time, "htrgate"
                    ).fits
                }
                if(!params.skiprfi) {
                    htrrfi_fits_path = "${params.publish_dir}/${params.label}/loadfits/htrrfi/${params.label}_htrrfi.fits"
                    if(new File(htrrfi_fits_path).exists()) {
                        htrrfi_fits = Channel.fromPath(htrrfi_fits_path)
                    }
                    else {
                        htrrfi_fits = corr_htrrfi(
                            "${params.label}_htrrfi", params.data_frb, params.ra_frb, 
                            params.dec_frb, opt_gate.htrrfi, opt_gate.polyco, 
                            binconfig.int_time, "htrrfi"
                        ).fits
                    }
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
