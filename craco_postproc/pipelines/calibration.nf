localise_dir = "$baseDir/../localise/"
beamform_dir = "$baseDir/../beamform/"

process determine_flux_cal_solns {
    publishDir "${params.publish_dir}/${params.label}/fluxcal", mode: "copy"

    input:
        path cal_fits
        val flagfile
        val target
        val cpasspoly

    output:
        path "calibration_noxpol_${target}.tar.gz", emit: solns
        path "*_calibrated_uv.ms", emit: ms
        path "*ps", emit: plots

    script:
        """
        if [ "$flagfile" == "" ]; then
            echo "You now need to write the flagfile for ${cal_fits} and provide it with --fluxflagfile!"
            exit 2
        fi        

        args="--calibrateonly"
        args="\$args -c $cal_fits"
        args="\$args --uvsrt"
        args="\$args -u 51"
        args="\$args --src=$target"
        args="\$args --cpasspoly=$cpasspoly"
        args="\$args -f 15"
        args="\$args --flagfile=$flagfile"

        $localise_dir/calibrateFRB.py \$args
        """
}

process apply_flux_cal_solns_finder {
    publishDir "${params.publish_dir}/${params.label}/finder", mode: "copy"

    input:
        path target_fits
        path cal_solns
        val target
        val cpasspoly

    output:
        // path "*.image", emit: all_images
        // path "f*.fits", emit: all_fits_images
        // path "*_calibrated_uv.ms", emit: all_ms
        path "*.jmfit", emit: all_jmfits
        // path "*.reg", emit: all_regs
        // path "${params.label}.image", emit: peak_image
        path "${params.label}.fits", emit: peak_fits_image
        path "${params.label}.jmfit", emit: peak_jmfit
        path "${params.label}.reg", emit: peak_reg
        path "${params.label}_calibrated_uv.ms", emit: peak_ms

    script:
        """
        tar -xzvf $cal_solns
        for b in `seq 0 19`; do
            bin="\$(printf "%02d" \$b)"

            args="--targetonly"
            args="\$args -t fbin\${bin}_norfi.fits"
            args="\$args -r 3"
            args="\$args --cpasspoly=$cpasspoly"
            args="\$args -i"
            args="\$args -j"
            args="\$args --cleanmfs"
            args="\$args --pols=I"
            args="\$args --imagename=finderbin\${bin}"
            args="\$args -a 16"
            args="\$args -u 500"
            args="\$args --skipplot"
            args="\$args --src=$target"
            args="\$args --nmaxsources=1"

            $localise_dir/calibrateFRB.py \$args

            for f in `ls finderbin\${bin}*jmfit`; do
                echo \$f
                $localise_dir/get_region_str.py \$f FRB >> finderbin\${bin}_sources.reg
            done
        done

        # Remove empty .jmfit and .reg files
        find . -type f -empty -print -delete

        # parse jmfits for S/N then find index of maximum
        SNs=`grep --no-filename "S/N" *jmfit | tr "S/N:" " "`
        peak_jmfit=\$($localise_dir/argmax.sh "\$(echo \$SNs)" "\$(ls *jmfit)")
        peak="\${peak_jmfit%.*}"
        peakbin=\${peak:9:2}

        echo "\$peak determined to be peak bin"
        cp \$peak_jmfit ${params.label}.jmfit
        # cp -r \${peak}.image ${params.label}.image
        cp \${peak}.fits ${params.label}.fits
        cp \${peak}_sources.reg ${params.label}.reg
        cp -r fbin\${peakbin}_norfi_calibrated_uv.ms ${params.label}_calibrated_uv.ms
        """    
}

process apply_flux_cal_solns_field {
    publishDir "${params.publish_dir}/${params.label}/field", mode: "copy"

    input:
        path target_fits
        path cal_solns
        val flagfile
        val target
        val cpasspoly
        val dummy   // so we can force this to wait for finder to finish

    output:
        // path "*.image", emit: image
        path "f*.fits", emit: fitsimage
        path "*_calibrated_uv.ms", emit: ms
        path "*jmfit", emit: jmfit
        path "*.reg", emit: regions

    script:
        """
        if [ "$flagfile" == "" ]; then
            echo "You now need to write the flagfile for ${target_fits} and provide it with --fieldflagfile!"
            exit 2
        fi    
        tar -xzvf $cal_solns

        args="--targetonly"
        args="\$args -t $target_fits"
        args="\$args -r 3"
        args="\$args --cpasspoly=$cpasspoly"
        args="\$args -i"
        args="\$args -j"
        args="\$args --cleanmfs"
        args="\$args --pols=I"
        args="\$args --imagename=field"
        args="\$args -a 16"
        args="\$args -u 500"
        args="\$args --skipplot"
        args="\$args --imagesize=3000"
        args="\$args --pixelsize=4"
        args="\$args --src=$target"
        args="\$args --tarflagfile=$flagfile"

        $localise_dir/calibrateFRB.py \$args
        i=1
        for f in `ls *jmfit`; do
            echo \$f
            $localise_dir/get_region_str.py \$f \$i >> sources.reg
            i=\$((i+1))
        done
        """    
}

process apply_flux_cal_solns_polcal {
    publishDir "${params.publish_dir}/${params.label}/polcal", mode: "copy"

    input:
        path target_fits
        path cal_solns
        val flagfile
        val target
        val cpasspoly

    output:
        // path "*.image", emit: image
        path "p*.fits", emit: fitsimage
        path "*_calibrated_uv.ms", emit: ms
        path "*jmfit", emit: jmfit
        path "*.reg", emit: regions

    script:
        """
        if [ "$flagfile" == "" ]; then
            echo "You now need to write the flagfile for ${target_fits} and provide it with --polcalflagfile!"
            exit 2
        fi    
        tar -xzvf $cal_solns

        args="--targetonly"
        args="\$args -t $target_fits"
        args="\$args -r 3"
        args="\$args --cpasspoly=$cpasspoly"
        args="\$args -i"
        args="\$args -j"
        args="\$args --cleanmfs"
        args="\$args --pols=I"
        args="\$args --imagename=polcal"
        args="\$args -a 16"
        args="\$args -u 500"
        args="\$args --skipplot"
        args="\$args --src=$target"
        args="\$args --tarflagfile=$flagfile"
        args="\$args --nmaxsources=1"

        $localise_dir/calibrateFRB.py \$args
        i=1
        for f in `ls *jmfit`; do
            echo \$f
            $localise_dir/get_region_str.py \$f \$i >> sources.reg
            i=\$((i+1))
        done
        """    
}

process determine_pol_cal_solns {
    publishDir "${params.publish_dir}/${params.label}/polcal", mode: "copy"
    
    input:
        path I
        path Q
        path U
        path V

    output:
        path "${params.label}_polcal.dat", emit: polcal_solns
        path "*.png", emit: plots
    
    script:
        """
        args="-i $I -q $Q -u $U -v $V"
        args="\$args -p $params.period_polcal"
        args="\$args -f $params.centre_freq_polcal"
        args="\$args -b 336"
        args="\$args -l ${params.label}_polcal"
        args="\$args -o ${params.label}_polcal.dat"
        args="\$args --reduce_df 1"
        args="\$args --plot"
        args="\$args --plotdir ."

        $beamform_dir/polcal.py \$args
        """
}
