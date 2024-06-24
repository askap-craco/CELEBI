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