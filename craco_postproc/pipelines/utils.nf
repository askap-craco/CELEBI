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
