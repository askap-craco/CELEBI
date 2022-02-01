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

process get_num_ant {
    input:
        val data

    output:
        stdout num_ants
    
    script:
        """
        ls -d $data/ak* | wc -w
        """
}
