/*
    Processes for calibration and imaging
*/

localise_dir = "$projectDir/../localise"
beamform_dir = "$projectDir/../beamform"

params.finderimagesize = 1024
params.finderpixelsize = 1
params.fieldimagesize = 5000
// new - adding in field pixel size
// params.fieldpixelsize = 1.5
if ( params.centre_freq_frb < 1000 ) {
	params.fieldpixelsize = 2.5
}
else if ( params.centre_freq_frb < 1200 ) {
	params.fieldpixelsize = 2
}
else {
	params.fieldpixelsize = 1.5
}

params.polcalimagesize = 128
params.minbeamfrac = 0.05
params.refant = 3   // reference antenna - index corresponds to ak name
                    // i.e. refant = 3 corresponds to ak03

params.nfieldsources = 50   // number of field sources to try and find
params.cpasspoly = 5
params.out_dir = "${params.publish_dir}/${params.label}"




process determine_flux_cal_solns {
    /*
        Determine flux calibration solutions

        Input
            cal_fits: path
                Flux calibrator visibilities in a FITS file
            flagfile: val
                Absolute path to AIPS flag file for flux calibrator
            fcm: path
                fcm file to update, unless already updated
            
        Output
            solns: path
                Tarball containing calibration solutions
            ms: path
                Calibrated flux calibrator measurement set
            plot: path
                AIPS postscript plots of calibration solutions
    */
    publishDir "${params.out_dir}/fluxcal", mode: "copy"

    input:
        path cal_fits
        val flagfile
        path fcm

    output:
        path "calibration_noxpol_${params.target}.tar.gz", emit: solns
        path "*_calibrated_uv.ms", emit: ms
        path "*ps", emit: plots
        path "fcm_delayfix.txt", emit: fcm_delayfix

    script:
        """
        ml apptainer
        set -a
        set -o allexport
        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"

        args="--calibrateonly"
        args="\$args -c $cal_fits"
        args="\$args --uvsrt"
        aipsid="\$((RANDOM%8192))"
        args="\$args -u \$aipsid"
        args="\$args --src=$params.target"
        args="\$args --cpasspoly=$params.cpasspoly"
        args="\$args -f 15"
        args="\$args --refant=$params.refant"
        if [ "$flagfile" != "" ]; then
            args="\$args --flagfile=$flagfile"
        fi
        # update fcm if not aready updated
        if [ "$fcm" != "fcm_delayfix.txt" ]; then
            cp $fcm fcm_delayfix.txt
            args="\$args --updatefcmfile=fcm_delayfix.txt"
        else
            touch fcm_delayfix.txt
        fi
        if [ "$params.nopossm" == "true" ]; then
            args="\$args --skipplot"
        fi

        apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ParselTongue $localise_dir/calibrateFRB.py \$args'
        rm -rf \$aips_dir 
        """
    
    stub:
        """
        touch calibration_noxpol_${params.target}.tar.gz
        touch stub_calibrated_uv.ms
        touch stub.ps
        touch fcm_delayfix.txt
        """
}

process image_finder {
    /*
        For a single finder bin:
            - Flux calibrate
            - Image
            - Find and fit a single source

        Input
            target_fits: path
                Finder bin visibilities in a FITS file
            cal_solns: path
                Tarball containing calibration solutions
        
        Output
            jmfit: path
                JMFIT outputs for this finder bin
            fits_image: path
                FITS-format image of this bin
            reg: path
                DS9 region file of source fit in this bin
            ms: path
                Calibrated bin measurement set
    */
    publishDir "${params.out_dir}/finder", mode: "copy"
    maxForks 1

    input:
        each path(target_fits)
        path cal_solns

    output:
        path "fbin*.jmfit", emit: jmfit
        path "fbin*.fits", emit: fits_image
        path "fbin*.reg", emit: reg
        path "*_calibrated_uv.ms", emit: ms

    script:
        """
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_parseltongue3
        #fi
        ml apptainer
        set -a
        set -o allexport
        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"

        tar -xzvf $cal_solns
        target_fits=$target_fits
        bin=\${target_fits:9:2}

        args="--targetonly"
        args="\$args -t $target_fits"
        args="\$args -r $params.refant"
        args="\$args -i"
        args="\$args -j"
        args="\$args --cleanmfs"
        args="\$args --pols=I"
        args="\$args --imagename=fbin\${bin}"
        args="\$args --imagesize=$params.finderimagesize"
        args="\$args --pixelsize=$params.finderpixelsize"
        args="\$args -a 16"
        aipsid="\$((RANDOM%8192))"
        args="\$args -u \$aipsid"
        args="\$args --skipplot"
        args="\$args --src=$params.target"
        args="\$args --nmaxsources=1"
        args="\$args --findsourcescript=$localise_dir/get_pixels_from_field.py"
        args="\$args --findsourcescript2=$localise_dir/get_pixels_from_field2.py"
        args="\$args --refant=$params.refant"

        if [ "$params.finderflagfile" != "" ] && [ "$params.finderflagfile" != "null" ]; then
            args="\$args --tarflagfile=$params.finderflagfile"
        fi

        apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ParselTongue $localise_dir/calibrateFRB.py \$args' 

        for f in `ls fbin\${bin}*jmfit`; do
            echo \$f
            apptainer exec $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/get_region_str.py \$f FRB >> fbin\${bin}_sources.reg'
        done
        rm -rf \$aips_dir
        """
        
    stub:
        """
        target_fits=$target_fits
        bin=\${target_fits:9:2}
        touch fbin\$bin.jmfit
        touch fbin\$bin.fits
        touch fbin\$bin.reg
        touch fbin\${bin}_calibrated_uv.ms
        """
}

process get_peak {
    /*
        Find the bin in which the S/N of the single fit source is highest,
        and return its related files. All outputs are renamed to [FRB name].*

        Input
            jmfit: path
                All JMFIT output files from finder bins
            fits_image: path
                All finder bin FITS images
            reg: path
                All finder bin DS9 region files
            ms: path
                All finder bin calibrated visibility measurement sets
        
        Output
            peak_jmfit: path
                JMFIT output file of peak bin
            peak_fits_image: path
                FITS image of peak bin
            peak_reg: path
                DS9 region file of peak bin
            peak_ms: path
                Calibrated visibility measurement set of peak bin    
    */
    publishDir "${params.out_dir}/finder", mode: "copy"

    input:
        path jmfit
        path fits_image
        path reg
        path ms
    
    output:
        path "${params.label}.jmfit", emit: peak_jmfit
        path "${params.label}.fits", emit: peak_fits_image
        path "${params.label}.reg", emit: peak_reg
        path "${params.label}_calibrated_uv.ms", emit: peak_ms

    script:
        """
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        # fi
        set -a
        set -o allexport

        # Remove empty .jmfit and .reg files
        find *jmfit -type f -empty -print -delete
        find *reg -type f -empty -print -delete

        #filter out jmfit files that have a large beam size
        BMINs=`grep --no-filename "Fit:" *jmfit | tr "x" " " | tr -d [:alpha:] | tr -d ':' | tr -d ';' | awk '{print \$1}'`
        BMAXs=`grep --no-filename "Fit:" *jmfit | tr "x" " " | tr -d [:alpha:] | tr -d ':' | tr -d ';' | awk '{print \$2}'`
        beamBMIN=`grep --no-filename "Fit:" *jmfit | tr "x" " " | tr -d [:alpha:] | tr -d ':' | tr -d ';' | awk '{print \$4}'`
        beamBMAX=`grep --no-filename "Fit:" *jmfit | tr "x" " " | tr -d [:alpha:] | tr -d ':' | tr -d ';' | awk '{print \$5}'`

        # largebeam_ind=\$(python3 $localise_dir/argBeamExceed.py "\$(echo \$BMINs)" "\$(echo \$BMAXs)" "\$(echo \$beamBMIN)" "\$(echo \$beamBMAX)" "\$(ls *jmfit)")

        ml apptainer
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/argBeamExceed.py "\$(echo \$BMINs)" "\$(echo \$BMAXs)" "\$(echo \$beamBMIN)" "\$(echo \$beamBMAX)" "\$(ls *jmfit)" >> largebeam_ind.txt'
        # echo "\$largebeam_ind are files to be removed"
        echo "files to be removed:"
        cat largebeam_ind.txt

        for file in \$(cat "largebeam_ind.txt")
        do
            mv \${file} \${file}REJECT
        done

        # parse jmfits for S/N then find index of maximum
        SNs=`grep --no-filename "S/N" *jmfit | tr "S/N:" " "`
        # peak_jmfit=\$(python3 $localise_dir/argmax.py "\$(echo \$SNs)" "\$(ls *jmfit)")

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/argmax.py "\$(echo \$SNs)" "\$(ls *jmfit)">>peak_jmfit.txt'
        peak_jmfit=\$(cat "peak_jmfit.txt")

        peak="\${peak_jmfit%.*}"
        peakbin=\${peak:4:2}

        echo "\$peak determined to be peak bin"
        cp \$peak_jmfit ${params.label}.jmfit
        cp \${peak}.fits ${params.label}.fits
        cp \${peak}_sources.reg ${params.label}.reg
        cp -r *bin\${peakbin}*calibrated_uv.ms ${params.label}_calibrated_uv.ms
        """    

    stub:
        """
        touch ${params.label}.jmfit
        touch ${params.label}.fits
        touch ${params.label}.reg
        touch ${params.label}_calibrated_uv.ms
        """
}

process image_field {
    /*
        Apply flux calibration to field visibilties and create field image,
        unless we have an already-made deep field image, then find and fit
        sources.

        Input
            target_fits: path
                Field visibilities in FITS file
            cal_solns: path
                Tarball containing calibration solutions
            flagfile: val
                Absolute path to AIPS flag file for field data
            dummy: val
                A dummy variable used to force field calibration to wait for
                finder calibration (otherwise calibrateFRB.py steps on itself)
        
        Output
            fitsimage: path, optional
                FITS format field image. Won't be output if using a deep field
                image
            ms: path, optional
                Calibrated field visibility measurement set. Won't be output if
                using a deep field image
            jmfit: path
                JMFIT output files for all sources found in image
            regions: path
                DS9 region file containing regions for all sources fit.
    */
    publishDir "${params.out_dir}/field", mode: "copy"

    input:
        path target_fits
        path cal_solns
        val flagfile
        val dummy

    output:
        path "f*.fits", emit: fitsimage, optional: true
        path "*_calibrated_uv.ms", emit: ms, optional: true
        path "*jmfit", emit: jmfit
        path "*.reg", emit: regions
        path "cutouts", emit: cutouts

    script:
        """
        tar -xzvf $cal_solns
        ml apptainer
        set -a
        set -o allexport
        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"

        args="--imagename=field"
        args="\$args -j"
        args="\$args -i"
        args="\$args --pols=I"
        aipsid="\$((RANDOM%8192))"
        args="\$args -u \$aipsid"
        args="\$args --imagesize=$params.fieldimagesize"
        args="\$args --findsourcescript=$localise_dir/get_pixels_from_field.py"
        args="\$args --findsourcescript2=$localise_dir/get_pixels_from_field2.py"
        args="\$args --nmaxsources=$params.nfieldsources"
        args="\$args --src=$params.target"
        args="\$args --refant=$params.refant"
        args="\$args --minbeamfrac=$params.minbeamfrac"

        # if we have an already-made field image, skip imaging
        if [ "$params.fieldimage" == "null" ]; then
            args="\$args --targetonly"
            args="\$args -t $target_fits"
            args="\$args --cleanmfs"
            args="\$args -a 16"
            args="\$args --skipplot"
            args="\$args --pixelsize=$params.fieldpixelsize"
            args="\$args --tarflagfile=$flagfile"
            if [ "$flagfile" != "" ]; then
                args="\$args --tarflagfile=$flagfile"
            fi
        else
            args="\$args --image=$params.fieldimage"
        fi

        apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ParselTongue $localise_dir/calibrateFRB.py \$args'
        i=1
        for f in `ls *jmfit`; do
            echo \$f
            apptainer exec $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/get_region_str.py \$f \$i >> sources.reg'
            i=\$((i+1))
        done
        apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ParselTongue $localise_dir/makefieldcutouts.py'
        rm -rf \$aips_dir
        """    
    
    stub:
        """
        touch f0.fits
        touch stub_calibrated_uv.ms
        touch stub.jmfit
        touch stub.reg
        """
}

process image_polcal {
    /*
        Apply flux calibration to and image polarisation calibrator
        visbilities, then fit a single source

        Input
            target_fits: path
                Polarisation calibrator visibilities in a FITS file
            cal_solns: path
                Tarball containing calibration solutions
            flagfile: val
                Absolute path to AIPS flag file for polarisation calibrator
        
        Output
            fitsimage: path
                FITS format image of polarisation calibrator
            ms: path
                Calibrated polarisation calibrator visibility measurement set
            jmfit: path
                JMFIT output for source fit
            regions: path
                DS9 region of source fit    
    */
    publishDir "${params.out_dir}/polcal", mode: "copy"

    input:
        path target_fits
        path cal_solns
        val flagfile

    output:
        path "p*.fits", emit: fitsimage
        path "*_calibrated_uv.ms", emit: ms
        path "*jmfit", emit: jmfit
        path "*.reg", emit: regions

    script:
        """
        tar -xzvf $cal_solns
        ml apptainer
        set -a
        set -o allexport
        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"

        args="--targetonly"
        args="\$args -t $target_fits"
        args="\$args -r 3"
        args="\$args -i"
        args="\$args -j"
        args="\$args --cleanmfs"
        args="\$args --pols=I"
        args="\$args --imagename=polcal"
        args="\$args --imagesize=$params.polcalimagesize"
        args="\$args -a 16"
        aipsid="\$((RANDOM%8192))"
        args="\$args -u \$aipsid"
        args="\$args --skipplot"
        args="\$args --src=$params.target"
        args="\$args --refant=$params.refant"
        if [ "$flagfile" != "" ]; then
            args="\$args --tarflagfile=$flagfile"
        fi
        args="\$args --nmaxsources=1"
        args="\$args --findsourcescript=$localise_dir/get_pixels_from_field.py"
        args="\$args --findsourcescript2=$localise_dir/get_pixels_from_field2.py"

        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_parseltongue3
        #
        # fi

        apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ParselTongue $localise_dir/calibrateFRB.py \$args'

        i=1
        for f in `ls *jmfit`; do
            echo \$f
            apptainer exec $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/get_region_str.py \$f \$i >> sources.reg'
            i=\$((i+1))
        done
        rm -rf \$aips_dir
        """

    stub:
        """
        touch p0.fits
        touch stub_calibrated_uv.ms
        touch stub.jmfit
        touch stub.reg
        """    
}

process image_htrgate {
    /*
        For a htrgate fits:
            - Flux calibrate
            - Image
            - Find and fit a single source

        Input
            target_fits: path
                htrgate bin visibilities in a FITS file
            cal_solns: path
                Tarball containing calibration solutions
        
        Output
            jmfit: path
                JMFIT output
            fits_image: path
                FITS-format image
            reg: path
                DS9 region file
            ms: path
                Calibrated measurement set
    */
    publishDir "${params.out_dir}/htrgate", mode: "copy"
    maxForks 1

    input:
        each path(target_fits)
        path cal_solns

    output:
        path "*.jmfit", emit: jmfit
        path "*.fits", emit: fits_image
        path "*.reg", emit: reg
        path "*_calibrated_uv.ms", emit: ms

    script:
        """
        # if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_parseltongue3
        # fi
        ml apptainer
        set -a
        set -o allexport
        aips_dir="/fred/oz313/tempaipsdirs/aips_dir_\$((RANDOM%8192))"
        cp -r /fred/oz313/aips-clean-datadirs \$aips_dir
        export APPTAINER_BINDPATH="/fred/oz313/:/fred/oz313/,\$aips_dir/DATA/:/usr/local/aips/DATA,\$aips_dir/DA00/:/usr/local/aips/DA00"

        tar -xzvf $cal_solns
        target_fits=$target_fits

        args="--targetonly"
        args="\$args -t $target_fits"
        args="\$args -r 3"
        args="\$args -i"
        args="\$args -j"
        args="\$args --cleanmfs"
        args="\$args --pols=I"
        args="\$args --imagename=fbin\${bin}"
        args="\$args --imagesize=$params.finderimagesize"
        args="\$args --pixelsize=$params.finderpixelsize"
        args="\$args -a 16"
        aipsid="\$((RANDOM%8192))"
        args="\$args -u \$aipsid"
        args="\$args --skipplot"
        args="\$args --src=$params.target"
        args="\$args --nmaxsources=1"
        args="\$args --findsourcescript=$localise_dir/get_pixels_from_field.py"
        args="\$args --findsourcescript2=$localise_dir/get_pixels_from_field2.py"
        args="\$args --refant=$params.refant"

        apptainer exec $params.container bash -c 'source /opt/setup_proc_container && ParselTongue $localise_dir/calibrateFRB.py \$args'

        for f in `ls *jmfit`; do
            echo \$f
            apptainer exec $params.container bash -c 'source /opt/setup_proc_container && python3 $localise_dir/get_region_str.py \$f FRB >> fbin\${bin}_sources.reg'

        done
        rm -rf \$aips_dir
        """
        
    stub:
        """
        touch htrgate.jmfit
        touch htrgate.fits
        touch htrgate.reg
        touch stub_calibrated_uv.ms
        """
}

process determine_pol_cal_solns {
    /*
        Determine polarisation calibration solutions.

        Input
            htr_data: path
                Stokes I, Q, U, V dynamic spectra of polarisation calibrator
                beamformed data as numpy files
        
        Output
            pol_cal_solns: path
                A file containing the delay (in ns) and phase offset solutions
                with errors
            plots: path
                A set of .png plots generated at various stages of polcal.py
                for troubleshooting/verifying solutions
    */
    publishDir "${params.out_dir}/polcal", mode: "copy"
    
    input:
        path htr_data

    output:
        path "${params.label}_polcal_solutions.txt", emit: pol_cal_solns
        path "*.png", emit: plots
    
    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
            #module load gcc/9.2.0
            #module load openmpi/4.0.2
            #module load python/3.7.4
            #module load numpy/1.18.2-python-3.7.4
            #module load matplotlib/3.2.1-python-3.7.4
            #module load scipy/1.6.0-python-3.7.4
            #module load astropy/4.0.1-python-3.7.4
        #fi
        ml apptainer
        set -a
        set -o allexport

        args="-i ${params.label}_polcal_I_dynspec_${params.dm_polcal}.npy"
        args="\$args -q ${params.label}_polcal_Q_dynspec_${params.dm_polcal}.npy"
        args="\$args -u ${params.label}_polcal_U_dynspec_${params.dm_polcal}.npy"
        args="\$args -v ${params.label}_polcal_V_dynspec_${params.dm_polcal}.npy"        
        args="\$args --l_model $params.polcal_l_model"
        args="\$args --v_model $params.polcal_v_model"
        args="\$args --priors $params.polcal_priors"

        args="\$args --peak_w $params.polcal_peak_w"
        args="\$args --rms_w $params.polcal_rms_w"
        args="\$args --tN $params.polcal_tN"
        args="\$args --fN $params.polcal_fN"
        args="\$args --RFIguard $params.polcal_guard"
        
        args="\$args --pa0 $params.polcal_pa0"
        args="\$args --f0 $params.polcal_f0"
        args="\$args --cfreq $params.centre_freq_polcal"
        args="\$args --bw $params.bw"

        args="\$args --cpus $params.polcal_cpus"
        args="\$args --live $params.polcal_live"

        if [ '$params.polcal_ellipse' == 'true' ]; then
            args="\$args --elipse"
        fi

        args="\$args --ofile ${params.label}_polcal_solutions.txt"

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/polcal.py \$args'
        
        cp polcal_sampler/polcal_corner.png .
        """
    
    stub:
        """
        touch ${params.label}_polcal_solutions.txt
        touch stub.png
        """
}


process apply_pol_cal_solns {
    /*
        
        Apply polcal solutions to FRB data
        
        Input:
            htr_path: path
                path to X and Y polarisation data for FRB
            polcal_solns: path
                full file path to polcal solutions

        Output:
            New X and Y data products with full polcal solutions applied


    */

    publishDir "${params.out_dir}/htr", mode: "copy"

    input: 
        val label
        path pol_time_series
        path pol_cal_solns
        val cfreq
        val dm

    output:
        path "*_calib_*.npy", emit: calib_data

    script:
        """
        ml apptainer
        set -a
        set -o allexport


        args="-x ${label}_X_t_${dm}.npy"
        args="\$args -y ${label}_Y_t_${dm}.npy"
        args="\$args --soln $pol_cal_solns"
        args="\$args --cfreq $cfreq"
        args="\$args --bw $params.bw"

        args="\$args --xout ${label}_calib_X_t_${dm}.npy"
        args="\$args --yout ${label}_calib_Y_t_${dm}.npy"

        if [ '$params.polcal_fast' == 'true' ]; then
            args="\$args --fast"
        fi

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $beamform_dir/apply_polcal.py \$args'

        """
    
    stub:
        """
        touch stub_calib_.npy
        """

}
