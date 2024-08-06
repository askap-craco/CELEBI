//
//	Processes for automated flagging
//

flagging_dir = "$baseDir/../flagging/"

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
        set -xu
	    badchanfile = ${flagging_dir}+"badchannels_askap_"+askapband+"_"+src+".txt"

        export ANKDIR=$params.ankdir

        python3 ${flagging_dir}/doflag.py ${infitsfile} outfitsfile.fits ${badchanfile} proper logfile.txt bad_ant_file.txt
        """
    
    stub:
        """
        touch test_flagging_proper.txt
        """
}






























