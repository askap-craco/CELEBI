//
//	Processes for automated flagging
//

flagging_dir = "$projectDir/../flagging"

//  Before calibration:
//
//	Run the process flag_proper on flux cal
//                  Inputs should be ***fluxcal.fits, bad_channel_file for cal
//		    Outputs are output_fitsfile, a log file and a list of bad antennas
//
//	Run the process flag_proper on pol cal
//                  Inputs should be ***polcal.fits, bad_channel_file for cal
//		    Outputs are output_fitsfile, a log file and a list of bad antennas
//
//	Run the process flag_initial on field
//                  Inputs should be ***field.fits, bad_channel_file for field
//		    Outputs are output_fitsfile, a log file
//
//  After calibration:
//
//	Run the process flag_proper on field
//                  Inputs should be ***field.fits, bad_channel_file for field
//		    Outputs are output_fitsfile, a log file and a list of bad antennas
//
//	List of bad channels for different bands are in the bad_channel_files
//	They are named as 'badchannels_askap_(low/mid/high)_(cal/field).txt'

process flag_proper {
    /*
        Flag visibilities 'properly'

        Input
            infitsfile: path    		(Input FITS file)     
	    askapband: val			(ASKAP band - low / mid / high)
	    src: val                            (cal / field)
			
        Output
            outfitsfile: path    		(Output FITS file)
            logfile: path               (Text file with logs/ may be useful if things go wrong)
            bad_ant_file: path          (List of bad antennas)
     
    */
    
    input:
        path infitsfile
        val askapband
	val src
    
    output:
        path 'outfitsfile.fits'
        path 'logfile.txt'
        path 'bad_ant_file.txt'
    
    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        #fi   
        
	badchanfile = ${flagging_dir}+"badchannels_askap_"+askapband+"_"+src+".txt"
	
	module load gcc/12.2.0 && module load gsl/2.7 && module load python/3.10.8 && module load numpy/1.24.2-scipy-bundle-2023.02 && module load matplotlib/3.7.0 && python3 ${flagging_dir}/doflag.py ${infitsfile} outfitsfile.fits ${badchanfile} proper logfile.txt bad_ant_file.txt
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
	    askapband: val			(ASKAP band - low / mid / high)
	    src: val                            (cal / field)
			
        Output
            outfitsfile: path    		(Output FITS file)
            logfile: path               (Text file with logs/ may be useful if things go wrong)
     
    */
    
    input:
        path infitsfile
        val askapband
	val src
    
    output:
        path 'outfitsfile.fits'
        path 'logfile.txt'
    
    script:
        """
        #if [ "$params.ozstar" == "true" ]; then
        #    . $launchDir/../setup_proc
        #fi   
        
	badchanfile = ${flagging_dir}+"badchannels_askap_"+askapband+"_"+src+".txt"
	
	module load gcc/12.2.0 && module load gsl/2.7 && module load python/3.10.8 && module load numpy/1.24.2-scipy-bundle-2023.02 && module load matplotlib/3.7.0 && python3 ${flagDir}/doflag.py ${infitsfile} outfitsfile.fits ${badchanfile} initital logfile.txt none
        """
    
    stub:
        """
        touch test_flagging_initial.txt
        """
}






























