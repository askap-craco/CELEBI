/*
    Processes for calibration and imaging
*/

localise_dir = "$projectDir/../localise"
beamform_dir = "$projectDir/../beamform"

params.finderimagesize = 1024
params.finderpixelsize = 1
params.fieldimagesize = 3000
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

    label 'python'
    label 'aips'

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
        source /opt/setup_proc_container
        set -x 

        aipsid="\$((RANDOM%8192))"

        args="--calibrateonly"
        args="\$args -c $cal_fits"
        args="\$args --uvsrt"
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
        # ensure casa is in the path
        export PATH=\$PATH:$params.casapath

        ParselTongue $localise_dir/calibrateFRB.py \$args
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

    label 'python'
    label 'aips'

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
        source /opt/setup_proc_container 

        aipsid="\$((RANDOM%8192))"

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

        ParselTongue $localise_dir/calibrateFRB.py \$args

        for f in `ls fbin\${bin}*jmfit`; do
            echo \$f
            source /opt/setup_proc_container && python3 $localise_dir/get_region_str.py \$f FRB >> fbin\${bin}_sources.reg
        done
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
    
    container "file://$params.container"

    label 'python'

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
        # Remove empty .jmfit and .reg files
        find *jmfit -type f -empty -print -delete
        find *reg -type f -empty -print -delete

        #filter out jmfit files that have a large beam size
        BMINs=`grep --no-filename "Fit:" *jmfit | tr "x" " " | tr -d [:alpha:] | tr -d ':' | tr -d ';' | awk '{print \$1}'`
        BMAXs=`grep --no-filename "Fit:" *jmfit | tr "x" " " | tr -d [:alpha:] | tr -d ':' | tr -d ';' | awk '{print \$2}'`
        beamBMIN=`grep --no-filename "Fit:" *jmfit | tr "x" " " | tr -d [:alpha:] | tr -d ':' | tr -d ';' | awk '{print \$4}'`
        beamBMAX=`grep --no-filename "Fit:" *jmfit | tr "x" " " | tr -d [:alpha:] | tr -d ':' | tr -d ';' | awk '{print \$5}'`

        source /opt/setup_proc_container 
        python3 $localise_dir/argBeamExceed.py \
                "\$(echo \$BMINs)" \
                "\$(echo \$BMAXs)" \
                "\$(echo \$beamBMIN)" \
                "\$(echo \$beamBMAX)" \
                "\$(ls *jmfit)" \
                >> largebeam_ind.txt

        echo "files to be removed:"
        cat largebeam_ind.txt

        for file in \$(cat "largebeam_ind.txt")
        do
            mv \${file} \${file}REJECT
        done

        # parse jmfits for S/N then find index of maximum
        SNs=`grep --no-filename "S/N" *jmfit | tr "S/N:" " "`

        python3 $localise_dir/argmax.py \
                "\$(echo \$SNs)" \
                "\$(ls *jmfit)" \
                >>peak_jmfit.txt

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

    label 'python'
    label 'aips'

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

    script:
        """
        source /opt/setup_proc_container

        aipsid="\$((RANDOM%8192))"

        tar -xzvf $cal_solns

        args="--imagename=field"
        args="\$args -j"
        args="\$args -i"
        args="\$args --pols=I"
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
            args="\$args -r 3"
            args="\$args --cleanmfs"
            args="\$args -a 16"
            args="\$args --skipplot"
            args="\$args --pixelsize=4"
            args="\$args --tarflagfile=$flagfile"
            if [ "$flagfile" != "" ]; then
                args="\$args --tarflagfile=$flagfile"
            fi
        else
            args="\$args --image=$params.fieldimage"
        fi

        ParselTongue $localise_dir/calibrateFRB.py \$args
        i=1
        for f in `ls *jmfit`; do
            echo \$f
            source /opt/setup_proc_container && python3 $localise_dir/get_region_str.py \$f \$i >> sources.reg
            i=\$((i+1))
        done
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
    label 'python'
    label 'aips'

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
        source /opt/setup_proc_container

        aipsid="\$((RANDOM%8192))"

        tar -xzvf $cal_solns
        
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


        ParselTongue $localise_dir/calibrateFRB.py \$args

        i=1
        for f in `ls *jmfit`; do
            echo \$f
            python3 $localise_dir/get_region_str.py \$f \$i >> sources.reg
            i=\$((i+1))
        done
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
    
    label 'python'
    label 'aips'

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
        source /opt/setup_proc_container

        aipsid="\$((RANDOM%8192))"

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
        args="\$args -u \$aipsid"
        args="\$args --skipplot"
        args="\$args --src=$params.target"
        args="\$args --nmaxsources=1"
        args="\$args --findsourcescript=$localise_dir/get_pixels_from_field.py"
        args="\$args --findsourcescript2=$localise_dir/get_pixels_from_field2.py"
        args="\$args --refant=$params.refant"

        ParselTongue $localise_dir/calibrateFRB.py \$args

        for f in `ls *jmfit`; do
            echo \$f
            python3 $localise_dir/get_region_str.py \$f FRB >> fbin\${bin}_sources.reg

        done
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

    label 'python'

    input:
        path htr_data

    output:
        path "${params.label}_polcal_solutions.txt", emit: pol_cal_solns
        path "*.png", emit: plots
    
    script:
        """
        source /opt/setup_proc_container 

        elipse=''
        if [ '$params.polcal_ellipse' == 'true' ]; then
            elipse="--elipse"
        fi


        python3 $beamform_dir/polcal.py \
                -i ${params.label}_polcal_I_dynspec_${params.dm_polcal}.npy \
                -q ${params.label}_polcal_Q_dynspec_${params.dm_polcal}.npy \
                -u ${params.label}_polcal_U_dynspec_${params.dm_polcal}.npy \
                -v ${params.label}_polcal_V_dynspec_${params.dm_polcal}.npy \
                --l_model $params.polcal_l_model \
                --v_model $params.polcal_v_model \
                --priors $params.polcal_priors \
                --peak_w $params.polcal_peak_w \
                --rms_w $params.polcal_rms_w \
                --tN $params.polcal_tN \
                --fN $params.polcal_fN \
                --RFIguard $params.polcal_guard \
                --pa0 $params.polcal_pa0 \
                --f0 $params.polcal_f0 \
                --cfreq $params.centre_freq_polcal \
                --bw $params.bw \
                --cpus $params.polcal_cpus \
                --live $params.polcal_live \
                --ofile ${params.label}_polcal_solutions.txt \
                \$elipse 
        
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
    
    label 'python'

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
        source /opt/setup_proc_container

        fast=' '
        if [ '$params.polcal_fast' == 'true' ]; then
            fast="--fast"
        fi

        python3 $beamform_dir/apply_polcal.py \
                -x ${label}_X_t_${dm}.npy \
                -y ${label}_Y_t_${dm}.npy \
                --soln $pol_cal_solns \
                --cfreq $cfreq \
                --bw $params.bw \
                --xout ${label}_calib_X_t_${dm}.npy \
                --yout ${label}_calib_Y_t_${dm}.npy \
                \$fast
        """
    
    stub:
        """
        touch stub_calib_.npy
        """

}
