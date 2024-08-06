/*
	Processes for automated flagging
	List of bad channels for different bands are in the bad_channel_files
	They are named as 'badchannels_askap_(low/mid/high)_(cal/field).txt'
*/

flagging_dir = "$projectDir/../flagging/"
params.out_dir = "${params.publish_dir}/${params.label}"

process flag_proper {
    /*
        Flag visibilities 'properly'

        Input
            infitsfile: path
                Input FITS file
	        outfitsfilepub: val
                Output FITS file name
	        src: val
                Mode, one of "cal" or "field"
			
        Output
	        outfile: path
                Flagged FITS file
            logfile: path
                Text file with logs
            bad_ant_file: path
                List of bad antennas
     
    */  
    label 'conda'
    
    input:
        path 'infitsfile.fits'
        val outfitsfilepub
        val src
    
    output:
	    path "outfitsfile.fits", emit: outfile
        path "logfile.txt", emit: logfile
        path "bad_ant_file.txt", emit: badant

    script:
        """	
        askap_band="low"

        if [ \$(echo "$params.centre_freq_frb > 1500.0" |bc -l) -gt 0 ] 
        then
            askap_band="high"
        elif [ \$(echo "$params.centre_freq_frb > 1000.0" |bc -l) -gt 0 ]
        then
        askap_band="mid"
        fi

        echo "FRB detected in ASKAP \$askap_band"

        badchanfile="${flagging_dir}badchannels_askap_\${askap_band}_${src}.txt"
        echo "Bad channel file \${badchanfile}"

        export ANKDIR=$params.ankdir
        /fred/oz313/anaconda/anaconda3/bin/conda run -n cracofunew python3 ${flagging_dir}doflag.py infitsfile.fits outfitsfile.fits \${badchanfile} proper logfile.txt bad_ant_file.txt

        cp outfitsfile.fits ${outfitsfilepub}
        """
    
    stub:
        """
        touch test_flagging_proper.txt
        """
}

process flag_initial {
    /*
        Flag known bad channels only 

        Input
            infitsfile: path
                Input FITS file
	        outfitsfile: val
                Output FITS file name
	        src: val
                Mode, one of "cal" or "field"
			
        Output
            outfile: path
                Flagged FITS file
            logfile: path
                Text file with logs
     
    */

    label 'conda'
    
    input:
        path infitsfile
        val outfitsfile
        val src
    
    output:
        path "${outfitsfile}.fits", emit: outfile
        path "logfile.txt", emit: logfile
    
    script:
        """
        
        askapband = low
            
        if [ $params.centre_freq_frb -gt 1500.0 ] 
        then
            askapband = high
        elif [ $params.centre_freq_frb -gt 1000.0 ] 
        then
            askapband = mid
        fi

        echo "FRB detected in ASKAP ${askapband}"

        badchanfile = ${flagging_dir}badchannels_askap_${askapband}_${src}.txt
        echo "Bad channel file ${badchanfile}"

        export ANKDIR=$params.ankdir
        /fred/oz313/anaconda/anaconda3/bin/conda run -n cracofunew python3 ${flagging_dir}doflag.py ${infitsfile} ${outfitsfile}.fits ${badchanfile} initital logfile.txt none
        """
    
    stub:
        """
        touch test_flagging_initial.txt
        """
}






























