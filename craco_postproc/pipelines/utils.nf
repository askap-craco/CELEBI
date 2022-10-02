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
