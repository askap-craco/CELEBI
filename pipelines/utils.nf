utils_dir = "$baseDir/../utils/"
params.out_dir = "${params.publish_dir}/${params.label}"

process create_empty_file {
    /*
        Creates an empty file of the desired filename

        Input
            filename: val
                Name of empty file to create
            
        Output
            file: path
                Empty file
    */
    cache 'lenient'

    input: 
        val filename

    output: 
        path "$filename"

    script:
        """
        touch $filename
        """
}

process do_filter_antenna {

    /*
        Check if each antenna has an X and Y pol with non-zero data, if not, filter out

        input:
            full list of antenna/pol dependant fine spectra with proper delays and calibrations,
            outputs of craftcor_tab.py
        output:
            list of antenna/pol files but only with antennas where both X and Y polarisations are present
            and non-zero

    */

    input:
        val label
        path ant

    output:
        path "*.npy"

    script:
        """
        ml apptainer
        set -a
        set -o allexport

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $utils_dir/filter_antenna.py'

        # save txt out file
        cp antenna_filtering.txt ${params.out_dir}/htr/info/${label}_antenna_filtering.txt
        """


}



workflow filter_antenna {

    /* 
    
        Take:
            list of file paths of all antenna/pol fine spectra
        Emit:
            filtered list of file paths where each antenna has non-zero X and Y pol products
    
    
    */

    take:
        label                   // label of process, frb or polcal
        unfiltered_files        // Full list of antena/pol filepaths, this is a tuple list with [pol, filepath]

    main:
        // filter out antenna data, this will give a list of files, both X and Y for each antenna that wasn't filtered
        do_filter_antenna(label, unfiltered_files.map {it[1]}.toList())

        // Need to output filterd files in same tuple format as inputs, we will use regexp
        filtered_Xfiles = Channel.of('X').combine(do_filter_antenna.out.flatten().filter(~/^.*(X_f_filtered.npy)$/))
        filtered_Yfiles = Channel.of('Y').combine(do_filter_antenna.out.flatten().filter(~/^.*(Y_f_filtered.npy)$/))

        filtered_files = filtered_Xfiles.concat(filtered_Yfiles)

    emit:
        filtered_ant = filtered_files
}





process compile_summary {

    publishDir "${params.out_dir}", mode: "copy"

    input:
        val _compile
        val _finalpos

    output:
        path "*.txt"

    script:
        """
        ml apptainer
        set -a
        set -o allexport

        args="-d ${params.out_dir}"
        args="\$args -l ${params.label}"
        args="\$args --cfreq ${params.centre_freq_frb}"
        args="\$args --bw ${params.bw}"

        ls ${params.out_dir}/*

        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 $utils_dir/compile_frb_summary.py \$args'
    
        """
}