process create_empty_file {
    input: 
        val filename

    output: 
        path "$filename"

    script:
        """
        touch $filename
        """
}

process get_num_ants {
    input:
        val data

    output:
        path "*.idx"
    
    script:
        """
        nants=\$(ls -d $data/ak* | wc -w)
        nantsm1=\$(( nants - 1 ))
        for i in \$(seq 0 \$nantsm1); do
            touch \${i}.idx
        done
        """
}
