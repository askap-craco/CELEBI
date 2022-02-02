process create_empty_file {
    input: 
        val filename

    output: 
        path $filename

    script:
        """
        touch $filename
        """
}

process get_num_ants {
    input:
        val data

    output:
        stdout
    
    script:
        """
        ls -d $data/ak* | wc -w
        """
}
