nextflow.enable.dsl=2
nextflow.enable.strict=true // be less generous

utils_dir    = "${projectDir}/../utils/"
beamform_dir = "${projectDir}/../beamform/"
localise_dir = "${projectDir}/../localise/"

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
        set -xu
        touch $filename
        """
}

process do_filter_antenna {

    /*
        Check if each antenna has an X and Y pol with non-zero data, if not, filter out

        input:
            full list of antenna/pol dependant fine spectra with proper delays and calibrations,
            outputs of craftcor_tab.py
        output:
            list of antenna/pol files but only with antennas where both X and Y polarisations are present
            and non-zero

    */	
	label 'celebi'
	
    input:
        val label
        path ant

    output:
        path "*.npy"

    script:
        """
        set -xu
        
        python3 $utils_dir/filter_antenna.py

        # save txt out file
        cp antenna_filtering.txt ${params.out_dir}/htr/info/${label}_antenna_filtering.txt
        """
}

workflow filter_antenna {

    /*     
        Take:
            list of file paths of all antenna/pol fine spectra
        Emit:
            filtered list of file paths where each antenna has non-zero X and Y pol products    
    */

    take:
        label                   // label of process, frb or polcal
        unfiltered_files        // Full list of antena/pol filepaths, this is a tuple list with [pol, filepath]

    main:
        // filter out antenna data, this will give a list of files, both X and Y for each antenna that wasn't filtered
        do_filter_antenna(label, unfiltered_files.map {it[1]}.toList())

        // Need to output filterd files in same tuple format as inputs, we will use regexp
        filtered_Xfiles = Channel.of('X').combine(do_filter_antenna.out.flatten().filter(~/^.*(X_f_filtered.npy)$/))
        filtered_Yfiles = Channel.of('Y').combine(do_filter_antenna.out.flatten().filter(~/^.*(Y_f_filtered.npy)$/))

        filtered_files = filtered_Xfiles.concat(filtered_Yfiles)

    emit:
        filtered_ant = filtered_files
}

process compile_summary {

    publishDir "${params.out_dir}", mode: "copy"
	
	label 'celebi'
	
    input:
        path outputs

    output:
        path "*.txt"

    script:
        """
        set -xu

        ls ${params.out_dir}/*

        python3 $utils_dir/compile_frb_summary.py \
        	-d ${params.out_dir} \
        	-l ${params.label} \
        	--cfreq ${params.centre_freq_frb} \
        	--bw ${params.bw} \
        	--dm ${params.dm_frb} \
        	-d ${params.out_dir}    
        """
}

process print_params {
    /*
        Prints the values of parameters used for this particular run    
    */
        
    publishDir "${params.out_dir}", mode: "copy"

    output:
        path "*.txt"

    script:
        """
        set -xu
        
        echo "\n****General****\n" > parameters.txt
        echo "label                   = $params.label" >> parameters.txt
        echo "celebi_container        = $params.celebi_container" >> parameters.txt
        echo "cracofunew_container    = $params.cracrofunew_container" >> parameters.txt
        echo "nbits                   = $params.nbits" >> parameters.txt
        echo "fcm                     = $params.fcm" >> parameters.txt	
        echo "snoopy                  = $params.snoopy" >> parameters.txt 
        
        echo "bw                      = $params.bw" >> parameters.txt 
        echo "uppersideband           = $params.uppersideband" >> parameters.txt 
        echo "hwfile                  = $params.hwfile" >> parameters.txt
        
        echo "cards                   = $params.cards" >> parameters.txt
        echo "fpgas                   = $params.fpgas" >> parameters.txt 

        echo "refant                  = $params.refant" >> parameters.txt
        echo "nants                   = $params.nants" >> parameters.txt
        echo "nants_fcal              = $params.nants" >> parameters.txt
        echo "nants_pcal              = $params.nants" >> parameters.txt
        echo "nants_frb               = $params.nants" >> parameters.txt

        echo "\n****Calibrators****\n" >> parameters.txt
        echo "data_fluxcal            = $params.data_fluxcal" >> parameters.txt
        echo "centre_freq_polcal      = $params.centre_freq_polcal" >> parameters.txt
        echo "data_polcal             = $params.data_polcal" >> parameters.txt
        echo "polcal_chanflag         = $params.polcal_chanflag" >> parameters.txt
        echo "polcal_name             = $params.polcal_name" >> parameters.txt
        echo "ra_polcal               = $params.ra_polcal" >> parameters.txt
        echo "dec_polcal              = $params.dec_polcal" >> parameters.txt
        echo "dm_polcal               = $params.dm_polcal" >> parameters.txt
        echo "polcal_F0               = $params.polcal_F0" >> parameters.txt
        echo "polcal_F1               = $params.polcal_F1" >> parameters.txt
        echo "polcal_MJD0             = $params.polcal_MJD0" >> parameters.txt
        echo "polcal_l_model          = $params.polcal_l_model" >> parameters.txt
        echo "polcal_v_model          = $params.polcal_v_model" >> parameters.txt
        echo "polcal_priors           = $params.polcal_priors" >> parameters.txt
        echo "polcal_pa0              = $params.polcal_pa0" >> parameters.txt
        echo "polcal_f0               = $params.polcal_f0" >> parameters.txt
        echo "polcal_peak_w           = $params.polcal_peak_w" >> parameters.txt
        echo "polcal_dynspec_sigma    = $params.polcal_dynspec_sigma" >> parameters.txt
        echo "polcal_baseline         = $params.polcal_baseline" >> parameters.txt
        echo "polcal_dynspec_tN       = $params.polcal_dynspec_tN" >> parameters.txt
        echo "polcal_dynspec_guard    = $params.polcal_dynspec_guard" >> parameters.txt
        echo "polcal_fN               = $params.polcal_fN" >> parameters.txt
        
        echo "\n****FRB****\n" >> parameters.txt
        echo "data_frb                = $params.data_frb" >> parameters.txt
        echo "ra_frb                  = $params.ra_frb"	>> parameters.txt
        echo "dec_frb                 = $params.dec_frb" >> parameters.txt
        echo "dm_frb                  = $params.dm_frb" >> parameters.txt
        echo "centre_freq_frb         = $params.centre_freq_frb" >> parameters.txt

        echo "\n****Localization****\n" >> parameters.txt        
        echo "bincinfig_gate          = $params.binconfig_gate " >> parameters.txt
        echo "polyco_gate             = $params.polyco_gate" >> parameters.txt
        echo "inttime_gate            = $params.inttime_gate" >> parameters.txt

        echo "usefield                = $params.usefield" >> parameters.txt                                                        
        echo "fieldimage              = $params.fieldimage" >> parameters.txt
        echo "numfinderbins           = $params.numfinderbins" >> parameters.txt
        echo "cenfinderbin            = $params.cenfinderbin" >> parameters.txt
        echo "searchms                = $params.searchms" >> parameters.txt
        echo "image_all_bins          = $params.image_all_bins" >> parameters.txt
        echo "fluxflagfile            = $params.fluxflagfile" >> parameters.txt
        echo "fieldflagfile           = $params.fieldflagfile" >> parameters.txt
        echo "polflagfile             = $params.polflagfile" >> parameters.txt
        echo "finderflagfile          = $params.finderflagfile" >> parameters.txt

        echo "\n****Beamforming****\n" >> parameters.txt
        
        echo "\n****Shrine****\n" >> parameters.txt
        
        echo "\n****MF Imaging****\n" >> parameters.txt

        """
}

