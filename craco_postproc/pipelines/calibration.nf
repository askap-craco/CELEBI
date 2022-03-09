localise_dir = "$baseDir/../localise/"

process determine_flux_cal_solns {
    publishDir "${params.publish_dir}/${params.label}/${task.process.replaceAll(':', '_')}", mode: "copy"

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

process apply_flux_cal_solns {
    publishDir "${params.publish_dir}/${params.label}/${task.process.replaceAll(':', '_')}", mode: "copy"

    input:
        path target_fits
        path cal_solns
        val flagfile
        val target
        val cpasspoly
        val dummy   // so we can force only one instance to go at a time

    output:
        path "*.image", emit: image
        path "f*.fits", emit: fitsimage
        path "*_calibrated_uv.ms", emit: ms
        path "*jmfit", emit: jmfit

    script:
        """
        if [ "$flagfile" == "" ]; then
            echo "You now need to write the flagfile for ${target_fits} and provide it with --polflagfile or --fieldflagfile!"
            exit 2
        fi    
        tar -xzvf $cal_solns

        if [ "$dummy" != "finder" ]; then
            args="--targetonly"
            args="\$args -t $target_fits"
            args="\$args -r 3"
            args="\$args --cpasspoly=$cpasspoly"
            args="\$args -i"
            args="\$args --dirtymfs"
            args="\$args --pols=I"
            args="\$args --imagename=field"
            args="\$args -a 16"
            args="\$args -u 500"
            args="\$args --skipplot"
            args="\$args --imagesize=2048"
            args="\$args --pixelsize=2"
            args="\$args --src=$target"

            if [ `wc -c $flagfile | awk '{print \$1}'` != 0 ]; then
                args="\$args --tarflagfile=$flagfile"
            fi

            $localise_dir/calibrateFRB.py \$args
            touch empty.stats
        else
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

                $localise_dir/calibrateFRB.py \$args
            done
        fi
        """    
}

process determine_pol_cal_solns {
    publishDir "${params.publish_dir}/${params.label}/${task.process.replaceAll(':', '_')}", mode: "copy"
    
    input:
        path htr_data

    output:
        path "polcal.dat"
    
    script:
        """
        touch polcal.dat
        """
}
