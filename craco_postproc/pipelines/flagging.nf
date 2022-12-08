//
//	Processes for automated flagging
//

flagging_dir = "$baseDir/../flagging/"

//  Before calibration:
//
//	Run the process flag_proper on flux cal
//                 
//
//	Run the process flag_proper on pol cal
//                  
//
//	Run the process flag_initial on field
//                 
//
//  After calibration:
//
//	Run the process flag_proper on field
//                
//
//	List of bad channels for different bands are in the bad_channel_files
//	They are named as 'badchannels_askap_(low/mid/high)_(cal/field).txt'

process flag_proper {
    /*
        Flag visibilities 'properly'

        Input
            infitsfile: path    		(Input FITS file)   
	    outfitsfile: val    		(Output FITS file)	       
	    src: val                            (cal / field)
			
        Output
	    outfile: path
            logfile: path              	        (Text file with logs/ may be useful if things go wrong)
            bad_ant_file: path                  (List of bad antennas)
     
    */
    
    input:
        path infitsfile
	val outfitsfile
	val src
    
    output:
	path ${outfitsfile}".fits", emit: outfile
        path "logfile.txt", emit: logfile
        path "bad_ant_file.txt", emit: badant
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi   
	
	if [ $params.centre_freq_frb -gt 1500.0 ]; then
	    askapband = "high"
        elif [ $params.centre_freq_frb -gt 1000.0 ]; then
            askapband = "mid"
        else
            askapband = "low"
        fi

        echo "FRB detected in ASKAP " ${askapband}

	badchanfile = ${flagging_dir}"badchannels_askap_"${askapband}"_"${src}".txt"
	echo "Bad channel file "$badchanfile

        python3 ${flagging_dir}/doflag.py ${infitsfile} ${outfitsfile}".fits" ${badchanfile} proper logfile.txt bad_ant_file.txt
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
            infitsfile: path    		(Input FITS file)   
	    outfitsfile: val    		(Output FITS file) 
	    src: val                            (cal / field)
			
        Output
	    outfile: path
            logfile: path                       (Text file with logs/ may be useful if things go wrong)
     
    */
    
    input:
        path infitsfile
	val outfitsfile
	val src
    
    output:
        path ${outfitsfile}".fits", emit: outfile
        path "logfile.txt", emit: logfile
    
    script:
        """
        if [ "$params.ozstar" == "true" ]; then
            . $launchDir/../setup_proc
        fi   

        if [ $params.centre_freq_frb -gt 1500.0 ]; then
	    askapband = "high"
        elif [ $params.centre_freq_frb -gt 1000.0 ]; then
            askapband = "mid"
        else
            askapband = "low"
        fi

        echo "FRB detected in ASKAP " ${askapband}

	badchanfile = ${flagging_dir}"badchannels_askap_"${askapband}"_"${src}".txt"
	echo "Bad channel file "$badchanfile

        python3 ${flagDir}/doflag.py ${infitsfile} ${outfitsfile}".fits" ${badchanfile} initital logfile.txt none
        """
    
    stub:
        """
        touch test_flagging_initial.txt
        """
}






























