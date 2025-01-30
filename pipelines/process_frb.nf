nextflow.enable.dsl=2
nextflow.enable.strict=true // be less generous

include { create_empty_file } from './utils'
include { correlate as corr_finder; correlate as corr_rfi;
    correlate as corr_gate; correlate as corr_field; 
    subtract_rfi as sub_rfi; subtract_rfi as sub_htrrfi; get_start_mjd as get_start_mjd } from './correlate'
include { image_finder; image_field; get_peak; image_htrgate } from './calibration'
include { find_offset; apply_offset; apply_offset as apply_offset_htr; 
    generate_binconfig } from './localise'
include { beamform as bform_frb; gen_dspec as frb_dspec; 
    dedisperse; ifft; generate_dynspecs } from './beamform'
include { flag_proper as flagdat } from './flagging'
include { shrine as smdm } from './shrine'
include { compile_summary } from './utils'
include { mf_image } from './mfimage'


utils_dir    = "${projectDir}/../utils/"
beamform_dir = "${projectDir}/../beamform/"
localise_dir = "${projectDir}/../localise/"


polarisations = Channel
    .fromList(params.pols)

antennas = Channel
    .of(0..params.nants_frb-1)

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
            fcm: path
                fcm file to use

        Output
            data: path
                1 ms time resolution dynamic spectrum
            time: path
                Time axis in MJD
    */

    label 'celebi'

    input:
        val label
        val data
        each pol
        each ant_idx
        path fcm

    output:
        path "${label}_ICS_${pol}*${ant_idx}.npy", emit: data
        path "t_mjd.npy", emit: time

    script:
        """
        export CRAFTCATDIR="."
        source /opt/setup_proc_container 
        set -xu

        startmjd=`python3 $localise_dir/get_start_mjd.py $data` 

        # Run processTimeStep.py with the --calconly flag to stop once calcfile is written
        python3 $localise_dir/processTimeStep.py \
                -t $data \
                --ra $params.ra_frb \
                -d $params.dec_frb" \
                -f $fcm \
                -b $params.nbits \
                --card 1 \
                -k \
                --name=${label}_ICS \
                -o . \
                --freqlabel c1_f0 \
                --dir=$localise_dir/../difx \
                --calconly \
                --startmjd \$startmjd

        mkdir delays    
        python3 $beamform_dir/craftcor_tab.py \
                -d $data \
                --parset $fcm \
                --calcfile c1_f0/craftfrb.im \
                -o ${label}_ICS \
                --ics \
                --cpus=8 \
                --pol=$pol \
                --an=$ant_idx
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
    publishDir "${params.publish_dir}/${params.label}/ics", mode: "copy"

    label 'celebi'

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
        source /opt/setup_proc_container
        set -xu
        
        python3 $localise_dir/sum_ics.py ${label}_ICS.npy ${label}_ICS_*.npy

        python3 $localise_dir/search_ics.py \
                --ds ${label}_ICS.npy \
                -s $snoopy \
                -t $t_mjd \
                -f $params.centre_freq_frb \
                --DMrange=$params.ICS_DMrange \
                --DMstep=$params.ICS_DMstep \
                -o ${label}.cand
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

    label 'celebi'

    output:
        env ra, emit: ra
        env dec, emit: dec

    script:
        """
        source /opt/setup_proc_container
        set -xu

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

        #export ant_pattern
        radec_beam=\$(python3 $localise_dir/get_beam_radec.py \$ra_beam_deg \$dec_beam_deg)
        echo \$radec_beam

        # radec_beam=`python3 $localise_dir/get_beam_radec.py \$ra_beam_deg \$dec_beam_deg`
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
            cand: path
                Refined candidate for FRB from ICS search
        
        Output:
            plot: path
                Plotted dynamic spectra across different time resolutions
            crops: path
                Directory containing cropped numpy files
    */
    publishDir "${params.publish_dir}/${params.label}/htr", mode: "copy"

    label 'celebi'

    input:
        val label
        path fnames_file
        path dynspecs
        val centre_freq
        val dm
        path cand
    
    output:
        path "*.png"
        path "crops", emit: crops
        path "crops/*.npy", emit: crop_us
        path "*IQUV*.png", emit: plot_file
    
    script:
        """
        source /opt/setup_proc_container
        set -xu
        
        start_time=`cat ${params.out_dir}/htr/info/bform_start_MJD.txt` 

        mkdir crops

        python3 $beamform_dir/plot.py \
                -s $fnames_file \
                -f $centre_freq \
                -l $label \
                -d $dm \
                -t \$start_time \
                -c $cand \
                --t_panels $params.plot_mosaic_t_list
        """
    
    stub:
        """
        touch stub.png
        mkdir crops
        touch crops/stub_50us_I.npy
        touch 50us_crop_start_s.txt
        """
}

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

    label 'celebi'

    input:
        path crops
        val dm

    output:
        env dmopt, emit: dm_opt
        path "*png"

    script:
        """
        source /opt/setup_proc_container
        set -xu

        python3 $beamform_dir/opt_DM.py \
                -x $crops/${params.label}_${dm}_X.npy \
                -y $crops/${params.label}_${dm}_Y.npy \
                -d $params.minDM \
                -D $params.maxDM \
                -s $params.DMstep \
                --DM0 $dm \
                --f0 $params.centre_freq_frb \
                --dt $params.opt_DM_dt 
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

    label 'celebi'

    input:
        path crop_50us
        path crop_start
    
    output:
        path "prof.txt", emit: prof

    script:
        """
        source /opt/setup_proc_container
        set -xu

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
        set -xu

        head -1 $polyco | awk '\$5="$dm"' > craftfrb.polyco
        head -2 $polyco | tail -1 | awk '\$6="1104.000"' >> craftfrb.polyco
        tail -1 $polyco >> craftfrb.polyco
        """
    
    stub:
        """
        cp old.polyco craftfrb.polyco
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
        fcm

    main:
    	
        if ( !params.skip_ics && params.nbits > 1 ) {
            coarse_ds = load_coarse_dynspec(params.label, params.data_frb, polarisations, antennas,fcm)
            refined_candidate_path = "${params.publish_dir}/${params.label}/ics/${params.label}.cand"            
            refine_candidate(params.label, coarse_ds.data.collect(), coarse_ds.time.first(), params.snoopy)
            refined_candidate = refine_candidate.out.cand           
        }
        else {
            refined_candidate = file(params.snoopy)
        }
    	
        field_fits_path = "${params.out_dir}/loadfits/field/${params.label}_field.fits"
        rfi_fits_path = "${params.out_dir}/loadfits/rfi/${params.label}_rfi.fits"
        finder_fits_path = "${params.out_dir}/loadfits/finder/finder*.fits"
        centre_bin_path = "${params.out_dir}/loadfits/finder/finderbin0${params.cenfinderbin}.fits"
        gate_fits_path = "${params.out_dir}/loadfits/gate/${params.label}_gate.fits"
        
        empty_file = create_empty_file("file")
    	
        if( params.localize || params.corrfrb ) {   
            
            binconfig = generate_binconfig(refined_candidate)      
            
            if(params.binconfig_gate != "") {
                binconfigpath = file(params.binconfig_gate).first()
                polycopath = file(params.polyco_gate).first()
                inttimepath = file(params.inttime_gate).first()
                // correlate gated FRB
                gate_fits = corr_gate("${params.label}_gate", params.data_frb, params.ra_frb, params.dec_frb,
                                       binconfigpath, polycopath, inttimepath, "gate", fcm).fits
            }
            else {
                // Correlate finder                
                (finder_fits, centre_bin_fits) = corr_finder(
                    "finder", params.data_frb, params.ra_frb, params.dec_frb, 
                    binconfig.finder, binconfig.polyco, binconfig.int_time, "finder", fcm
                )                

                // Correlate RFI (if not directly flagging finder)                
                if(!params.skiprfi) {
                    rfi_fits = corr_rfi(
                        "${params.label}_rfi", params.data_frb, params.ra_frb, 
                        params.dec_frb, binconfig.rfi, binconfig.polyco, binconfig.int_time, "rfi",
                        fcm
                    ).fits
                }
            }
        }
        else {
            finder_fits = file(finder_fits_path)
            centre_bin_fits = file(centre_bin_path)
            rfi_fits = file(rfi_fits_path)
            gate_fits = file(gate_fits_path)
        }
        
        if( params.localize || params.corrfld ) {
            
            if( !params.usefield ) {
                // Correlate field (if not using deep field image)            
                beam_centre = get_beam_centre()
                field_fits = corr_field(
                    "${params.label}_field", params.data_frb, beam_centre.ra, 
                    beam_centre.dec, empty_file, empty_file, empty_file, "field", fcm
                ).fits
            }
            else {
                field_fits = file(field_fits_path)
            }
        }
        else {
            field_fits = file(field_fits_path)
        }

        // Imaging field
        fld_srcs_path =  "${params.out_dir}/field/*.jmfit"
        if( params.localize || params.imgfld) {
            // Flagging
            if( !params.noflag && !params.usefield ) {
                field_fits_flagged = "${params.out_dir}/loadfits/field/${params.label}_field_f.fits"
                field_outfits = flagdat(field_fits,field_fits_flagged, "field").outfile           
                field_fits = field_outfits
            }

            field_sources = image_field(
                field_fits, flux_cal_solns, params.fieldflagfile, "NULL"
            ).jmfit        
        }        
        else {
            field_sources = file(fld_srcs_path)
        }

        frb_jmfit_path = "${params.out_dir}/finder/${params.label}.jmfit"
        // Imaging FRB        
        if( params.localize || params.imgfrb) {  
            
            if( params.imgfrb ) {
                binconfig = generate_binconfig(refined_candidate)
            }

            if(params.binconfig_gate != ""){
                gate_out = image_htrgate(gate_fits, flux_cal_solns)
                gate_jmfits = gate_out.jmfit
                gate_fits_images = gate_out.fits_image
                gate_regs = gate_out.reg
                gate_mss = gate_out.ms

                askap_frb_pos = get_peak(
                    gate_jmfits.collect(), gate_fits_images.collect(),
                    gate_regs.collect(), gate_mss.collect()
                ).peak_jmfit
            }
            else {
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
                bin_mss = bins_out.mstar

                askap_frb_pos = get_peak(
                    bin_jmfits.collect(), bin_fits_images.collect(), 
                    bin_regs.collect(), bin_mss.collect()
                ).peak_jmfit
            }
        }
        else {
            askap_frb_pos = file(frb_jmfit_path)
        }

        // Get FRB position        
        if( params.localize || params.getpos) {  
                
            offset_path = "${params.out_dir}/position/offset0.dat"
	        doffset_path = "${params.out_dir}/position/offsetfit.txt"
            frb_pos_path = "${params.out_dir}/position/${params.label}_final_position.txt"
                          
    		offres = find_offset(field_sources)
            offset = offres.offset
            doffset = offres.doffset                

        	finalres = apply_offset(offset, doffset, askap_frb_pos)
            final_position = finalres.final_position
        	// finalmap = finalres.hpmap
        }
            	
        final_position_path = "${params.out_dir}/finder/${params.label}_final_position.txt"
        final_position = file(final_position_path)



        // Beamforming
        xypath = "${params.out_dir}/htr/${params.label}_*_t_${params.dm_frb}.npy"
        if( params.gethtr || params.beamfrb ) {
            bform_frb(
                params.label, params.data_frb, askap_frb_pos, flux_cal_solns, 
                pol_cal_solns, params.dm_frb, params.centre_freq_frb,
                params.nants_frb, fcm, params.snoopy
            )
        }
        else {
            xy = file(xypath)
        }
        
        // Generate dynamic spectra
        if( params.gethtr || params.plotfrb ) {
            frb_dspec(
                params.label, xy, params.dm_frb, params.centre_freq_frb
            )                              
            plot(
                params.label, frb_dspec.out.dynspec_fnames, frb_dspec.out.htr_data,
                params.centre_freq_frb, params.dm_frb, refined_candidate
            )
        }
        
                
        if( params.shrine ) {
        	
        	idspath	= "${params.out_dir}/htr/crops/${params.label}_${params.dm_frb}_dsI_crop.npy"        	
        	idsdata	= file(idspath)
        	
        	smdm(idsdata,params.timresus)        	
        }             


        if( params.mfimage ) {

            // paths to required files
            ids_path = file("${params.out_dir}/htr/crops/${params.label}_${params.dm_frb}_dsI_crop.npy")
            binconfig = file("${params.out_dir}/binconfigs/craftfrb.finder.binconfig")
            polyco = file("${params.out_dir}/binconfigs/craftfrb.polyco")
            summary = file("${params.out_dir}/${params.label}_summary.txt")

            // do mf imaging
            mf_final_position = mf_image(ids_path, binconfig, polyco, 
                                    summary, flux_cal_solns, fcm).mf_final_position

        }


        // Code to compile output files into a single summary .txt file
        compile_out = Channel.empty()

        // If Imaging and Beamforming is being done
        if ( ( params.localize || params.getpos ) & (params.beamfrb || params.gethtr ) ) {
            params.do_compile_summary = true
            compile_out = compile_out.concat(finalres.final_position, plot.out.plot_file)


        } // else if only imaging is being done
        else if ( params.localize || params.getpos ) {
            params.do_compile_summary = true
            compile_out = compile_out.concat(finalres.final_position)

        } // else if only beamforming is being done
        else if ( params.plotfrb || params.gethtr ) {
            params.do_compile_summary = true
            compile_out = compile_out.concat(plot.out.plot_file)

        } // else if matched filter imaging is being done
        else if ( params.mfimage ) {
            params.do_compile_summary = true
            compile_out = compile_out.concat(mf_final_position)

        }
        else { // skip if none of the above
            params.do_compile_summary = false
        }

        if ( params.do_compile_summary ) {

            compile_summary(compile_out)
        }
}
