process create_empty_file {
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
