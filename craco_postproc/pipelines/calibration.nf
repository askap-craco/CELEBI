process determine_flux_cal_solns {
    input:
        path cal_fits
        path flagfile
        val target
        val cpasspoly

    output:
        path "calibration_noxpol_${target}.tar.gz"

    script:
        """
        args="--calibrateonly"
        args="\$args -c $cal_fits"
        args="\$args --uvsrt"
        args="\$args -u 51"
        args="\$args --src=$target"
        args="\$args --cpasspoly=$cpasspoly"
        args="\$args -f 15"
        args="\$args --flagfile=$flagfile"

        calibrateFRB.py \$args
        """
}

process apply_flux_cal_solns {
    input:
        path target_fits
        path cal_solns
        path flagfile
        val target
        val cpasspoly

    output:
        path "*.image"

    script:
        """
        tar -xzvf $cal_solns

        args="--targetonly"
        args="\$args -t $target_fits"
        args="\$args -r 3"
        args="\$args --cpasspoly=$cpasspoly"
        args="\$args -i"
        args="\$args --dirtymfs"
        args="\$args -a 16"
        args="\$args -u 500"
        args="\$args --skipplot"
        args="\$args --tarflagfile=$flagfile"
        args="\$args --src=$target"

        calibrateFRB.py \$args
        """    
}

process determine_pol_cal_solns {
    input:
        path htr_data

    output:
        path "polcal.dat"
}
