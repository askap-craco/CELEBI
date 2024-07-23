nextflow.enable.dsl=2

// Strictly FRB related params
params.dm_low = -5
params.dm_high = 5
params.dm_step = 0.1
//these are wide and coarse values, respectively
params.dm_count = 0
params.crop_dur = 100   // length of crop in ms
params.force_peak = 0   // force peak to a given timestep (ms). 0 to disable

params.timescale = 1		//	SMDM time resolution in us
params.bandwidth = 336      //  Bandwidth in MHz

//Debugging/diagnosis params
params.force_kc = 0
params.do_vary_kc = false
params.do_sn = false
params.do_uncertainty_min = false

//Output Params
params.saving = true

shrine_dir = "$projectDir/../shrine"
params.out_dir = "${params.publish_dir}/${params.label}"

process generate_profiles {
    /*
        Generates intensity profiles for structure maximisation.

        Input
            label: val
                FRB name
            dm: val
                Dispersion measure to which the input has been dedispersed in pc/cm3
            dm_low: val
                Lower end of desired delta DM range
            dm_high: val
                Upper end of desired delta DM range
            dm_step: val
                Delta DM step size in pc/cm3
            dm_count: val
                Number of samples to divide delta DM range into. Optional, overrides dm_step if used.
            timescale: val
                Time resolution to return in us
            centre_freq: val
                Central frequency of observation in MHz
            bandwidth: val
                Bandwidth of observation in MHz
            indspec: path
                Input intensity dynamic spectrum
            tresus: val
                Input time resolution in us

        Output
            DMdata: path
                Set of delta DMs corresponding to Idata
            Idata: path
                I(DM,t) profile for the FRB
            *.dat: path
                Dat versions of the above for backwards compatibility with earlier code.
            summary: path
                File summarising outputs of the script that may not appear in plots
    */    
    publishDir "${params.out_dir}/shrine", mode: "copy"

    input:
        val label
        val dm
        //val dm_range
        val dm_low
        val dm_high
        val dm_step
        val dm_count
        val timescale
        val centre_freq
        val bandwidth
        path indspec
        val tresus

    output:
        path("${label}_DMs.npy"), emit: DMdata
        path("${label}_I_${timescale}us.npy"), emit: Idata
        path("*.dat")
        path("${label}_profile_summaryfile.txt"), emit: summary

    script:
        """
        ml apptainer
        set -a
        set -o allexport
        
        args="-l ${label}"
        args="\$args -d ${dm}"
        args="\$args -L ${dm_low}"
        args="\$args -H ${dm_high}"
        args="\$args --dDM ${dm_step}"
        args="\$args --cDM ${dm_count}"
        args="\$args -t ${timescale}"
        args="\$args -f ${centre_freq}"
        args="\$args --crop_dur ${params.crop_dur}"
        args="\$args -I ${indspec}"
	    args="\$args -t0 ${tresus}"
        if [ "${params.force_peak}" != "0" ]; then
            args="\$args --force_peak ${params.force_peak}"
        fi
        args="\$args --bw ${params.bandwidth}"
		
		echo "python3 ${shrine_dir}/generate_profiles.py \$args > ${params.out_dir}/shrine/shrinelog.out"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 ${shrine_dir}/generate_profiles.py \$args > ${params.out_dir}/shrine/shrinelog.out'
        
        """

    stub:
        """
        touch ${label}_DMs.npy
        touch ${label}_I_${timescale}us.npy
        touch profiles.dat
        touch ${label}_profile_summaryfile.txt
        """

}


process qlookplot {
    /*
        Generates quick look plots for structure maximized FRB profile

        Input
            i_ds: path
            	Input Stokes I dynamic spectrum
            tresus: val
            	Input time resolution in us
            smresfile: path
            	File containing structure maximization results
        Output
            idsmdm: path
                Structure maximized Intensity plots
    */
    
    publishDir "${params.out_dir}/shrine", mode: "copy"
																
    input:
        path i_ds
        val tresus
        path smresfile

    output:
        path("${params.label}_idsmdm.png"), emit: idsmdm

    script:
        """
        ml apptainer
        set -a
        set -o allexport
        
        args="-l ${params.label}"
        args="\$args -t0 ${tresus}"
        args="\$args -t ${params.timescale}"
        args="\$args -f ${params.centre_freq_frb}"
        args="\$args --bw ${params.bandwidth}"
        args="\$args -I ${i_ds}"
        args="\$args -S ${smresfile}"        
			
        echo "python3 ${shrine_dir}/qlookplot.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 ${shrine_dir}/qlookplot.py \$args'
        
        """

    stub:
        """
        touch ${params.label}_idsmdm.png
        """
}


process maximise_structure {
    /*
        Determines the structure maximising delta DM and associated uncertainty

        Input
            label: val
                FRB name
            dm: val
                Dispersion measure to which the input has been dedispersed in pc/cm3
            timescale: val
                Time resolution to return in us
            saving: bool
                Whether or not to save plots.
            force_kc: val
                A value of kc to use in lieu of calculating it. Optional
            DMdata: path
                Set of delta DMs corresponding to Idata
            Idata: path
                I(DM,t) profile for the FRB
        
        Output:
            *.dat: path
                Files containing the structure parameter, uncertainty and relative uncertainty values. (all against DM)
            *.png: path
                Various plots, see README
            summary: path
                File summarising outputs of the script that may not appear in plots
            DM: env
                Structure maximising delta DM index
            resfile: path
            	File containing structure maximized DM
    */

    publishDir "${params.out_dir}/shrine", mode: "copy"

    input:
        val label
        val dm
        val timescale
        val saving
        val force_kc
        path DMdata
        path Idata

    output:
        path("*.dat")
        path("${label}*.png"), optional: true
        path("${label}_structure_summaryfile.txt"), emit: summary
        env DM, emit: DM, optional: true
        path("${label}_smdm_result.npy"), emit:resfile

    script:
        """
        ml apptainer
        set -a
        set -o allexport
        
        args="-l ${label}"
        args="\$args -d ${dm}"
        args="\$args -t ${timescale}"
        if [ "${saving}" == "true" ]; then
            args="\$args -s"   
        fi
        if [ "${force_kc}" != "0" ]; then
            args="\$args -kc ${force_kc}"   
        fi
				
		echo "python3 ${shrine_dir}/maximise_structure.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 ${shrine_dir}/maximise_structure.py \$args'
		
        if [ "${saving}" == "true" ]; then
            DM=\$( cat DM.txt )
        fi

        """

    stub:
        """
        touch ${label}_structure_summaryfile.txt
        touch structure.dat
        if [ "${saving}" == "true" ]; then
            DM=1
        fi
        """

}

process maximise_sn {
    /*
        Determines the signal to noise maximising delta DM

        Input
            label: val
                FRB name
            dm: val
                Dispersion measure to which the input has been dedispersed in pc/cm3
            timescale: val
                Time resolution to return in us
            saving: bool
                Whether or not to save plots.
            DMdata: path
                Set of delta DMs corresponding to Idata
            Idata: path
                I(DM,t) profile for the FRB
        
        Output:
            *.png: path
                Plots, see README
            summary: path
                File summarising outputs of the script that may not appear in plots
            ${label}_SN_DMs.dat: path
                Adjusted delta DM set with 0 as the signal to noise maximising DM.
            Ismooth: path
                Smoothed I(DM,t) profile to hand to structure_SN_comparison
            DM: env
                Signal to noise ratio maximising delta DM index
    */

    publishDir "${params.out_dir}/shrine", mode: "copy"

    input:
        val label
        val timescale
        val saving
        val force_kc
        path DMdata
        path Idata

    output:
        path("*.png"), optional: true
        path("${label}_SN_summaryfile.txt"), emit: summary
        path("${label}_SN_DMs.dat")
        path("${label}_I_smooth.npy"), emit: Ismooth, optional: true
        env DM, emit: DM, optional: true

    script:
        """
        ml apptainer
        set -a
        set -o allexport
        
        args="-l ${label}"
        args="\$args -t ${timescale}"
        if [ "${saving}" == "true" ]; then
            args="\$args -s"   
        fi
        if [ "${force_kc}" != "0" ]; then
            args="\$args -kc ${force_kc}"   
        fi
		        
        echo "python3 ${shrine_dir}/maximise_sn.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 ${shrine_dir}/maximise_sn.py \$args'

        if [ "${saving}" == "true" ]; then
            DM=\$( cat DM.txt )
        fi

        """

    stub:
        """
        touch ${label}_SN_summaryfile.txt
        touch ${label}_SN_DMs.dat
        ${label}_I_smooth.npy
        if [ "${saving}" == "true" ]; then
            DM=1
        fi
        """

}


process compare_structure_SN {
    /*
        Plots a comparison of the smoothed signal at the structure and S/N maximising delta DMs.

        Input
            label: val
                FRB name
            timescale: val
                Time resolution to return in us
            saving: bool
                Whether or not to save plots.
            DMdata: path
                Set of delta DMs corresponding to Idata
            Idata: path
                I(DM,t) profile for the FRB
            Ismooth: path
                Smoothed I(DM,t) profile for the FRB
            structure_DM: val
                Index of the structure maximising DM in DMdata
            SN_DM: val
                Index of the signal to noise ratio maximising DM in DMdata
        
        Output:
            *.png: path
                Plot, see README
    */

	publishDir "${params.out_dir}/shrine", mode: "copy"

    input:
        val label
        val timescale
        path DMdata
        path Idata
        path Ismooth
        val structure_DM
        val SN_DM

    output:
        path("*.png")

    script:
        """
        ml apptainer
        set -a
        set -o allexport
        
        args="-l ${label}"
        args="\$args -t ${timescale}"
        args="\$args -i_st ${structure_DM}"
        args="\$args -i_sn ${SN_DM}"
        args="\$args -d ${DMdata}"
        args="\$args -I ${Idata}"
        args="\$args -s ${Ismooth}"
		        
        echo "python3 ${shrine_dir}/structure_SN_comparison.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 ${shrine_dir}/structure_SN_comparison.py \$args'

        """

    stub:
        """
        touch compare.png
        """
}


process vary_kc {
    /*
        Generates plots showing the effect of different kc values on the structure parameter and delta DM

        Input
            label: val
                FRB name
            dm: val
                Dispersion measure to which the input has been dedispersed in pc/cm3
            timescale: val
                Time resolution to return in us
            saving: bool
                Whether or not to save plots.
            force_kc: val
                A value of kc to use in lieu of calculating it. Optional
            DMdata: path
                Set of delta DMs corresponding to Idata
            Idata: path
                I(DM,t) profile for the FRB
        
        Output:
            *.png: path
                Plot, see README
    */

    publishDir "${params.out_dir}/shrine", mode: "copy"

    input:
        val label
        val timescale
        val saving
        val force_kc
        path DMdata
        path Idata

    output:
        path("*.png")

    script:
        """
        ml apptainer
        set -a
        set -o allexport
        
        args="-l ${label}"
        args="\$args -t ${timescale}"
        if [ "${saving}" == "true" ]; then
            args="\$args -s"   
        fi
        if [ "${force_kc}" != "0" ]; then
            args="\$args -kc ${force_kc}"   
        fi
		        
        echo "python3 ${shrine_dir}/vary_kc.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 ${shrine_dir}/vary_kc.py \$args'

        """

    stub:
        """
        touch vary_kc.png
        """

}

process minimise_uncertainty {
    /*
        Determines the kc value which minimises uncertainty
        Generates a plot showing the effects of kc variation on delta DM and uncertainty.

        Input
            label: val
                FRB name
            dm: val
                Dispersion measure to which the input has been dedispersed in pc/cm3
            timescale: val
                Time resolution to return in us
            saving: bool
                Whether or not to save plots.
            force_kc: val
                A value of kc to use in lieu of calculating it. Optional
            DMdata: path
                Set of delta DMs corresponding to Idata
            Idata: path
                I(DM,t) profile for the FRB
        
        Output:
            *.png: path
                Plots, see README
            kc: env
                Calculated uncertainty minimising kc value.
            summary: path
                File summarising outputs of the script that may not appear in plots
    */

    publishDir "${params.out_dir}/shrine", mode: "copy"

    input:
        val label
        val timescale
        val saving
        val force_kc
        path DMdata
        path Idata

    output:
        path("*.png")
        env kc, emit: kc
        path("${label}_uncertainty_summaryfile.txt"), emit: summary

    script:
        """
        ml apptainer
        set -a
        set -o allexport
        
        args="-l ${label}"
        args="\$args -t ${timescale}"
        if [ "${saving}" == "true" ]; then
            args="\$args -s"   
        fi
        if [ "${force_kc}" != "0" ]; then
            args="\$args -kc ${force_kc}"   
        fi
		        
        echo "python3 ${shrine_dir}/minimise_uncertainty.py \$args"
        apptainer exec -B /fred/oz313/:/fred/oz313/ $params.container bash -c 'source /opt/setup_proc_container && python3 ${shrine_dir}/minimise_uncertainty.py \$args'

        kc=\$( cat kc.txt )
        """

    stub:
        """
        touch minimise_uncertainty.png
        kc=1
        touch ${label}_uncertainty_summaryfile.txt
        """

}

process make_summary{
    /*
        Generates a summary file containing all provided parameters.

        Input
            Technically, none. Does need parameters to exist.
        
        Output:
            summary: path
                File summarising various input parameters.
    */
    
    publishDir "${params.out_dir}/shrine", mode: "copy"
    
    output:
        path("full_summary.txt"), emit: summary

    script:
        """        
        echo "//User inputs:" > full_summary.txt
        echo "params.label = ${params.label}" >> full_summary.txt
        echo "params.dm_low = ${params.dm_low}" >> full_summary.txt
        echo "params.dm_high = ${params.dm_high}" >> full_summary.txt

        if [ "${params.dm_count}" == "0" ]; then
            echo "params.dm_step = ${params.dm_step}" >> full_summary.txt
            echo "//params.dm_count = ${params.dm_count}\t//default value so dm_step was used instead" >> full_summary.txt
        else
            echo "//params.dm_step = ${params.dm_step}\t//dm_count was supplied so this value wasn't used" >> full_summary.txt
            echo "params.dm_count = ${params.dm_count}" >> full_summary.txt
        fi
        echo "params.timescale = ${params.timescale}" >> full_summary.txt
        echo "params.crop_dur = ${params.crop_dur}" >> full_summary.txt
        echo "params.bandwidth = ${params.bandwidth}" >> full_summary.txt

        if [ "${params.force_kc}" == "0" ]; then
            echo "params.force_kc = ${params.force_kc} \t//default value so no forced kc value was used." >> full_summary.txt
        else
            echo "params.force_kc = ${params.force_kc}" >> full_summary.txt
        fi
        echo "params.do_vary_kc = ${params.do_vary_kc}" >> full_summary.txt
        echo "params.do_sn = ${params.do_sn}" >> full_summary.txt
        echo "params.do_uncertainty_min = ${params.do_uncertainty_min}" >> full_summary.txt
        echo "params.saving = ${params.saving}" >> full_summary.txt

        echo "//params taken from nextflow.config:" >> full_summary.txt
        echo "params.configs = ${params.configs}" >> full_summary.txt
        echo "params.data = ${params.data}" >> full_summary.txt

        echo "//params taken from external config file:" >> full_summary.txt
        echo "//looked for config file at ${params.configs}/${params.label}.config" >> full_summary.txt
        echo "params.dm_frb = ${params.dm_frb}" >> full_summary.txt
        echo "params.centre_freq_frb = ${params.centre_freq_frb}" >> full_summary.txt

        echo "" >> full_summary.txt
        """
}

process cat_summaries{
    /*
        Appends a text file or set of text files to another text file..

        Input
            main: path
                The file onto which others should be appended.
            summary: path
                The file(s) to append.
        
        Output:
            full_summary.txt: path
                File with other files' contents appended
    */
    
    publishDir "${params.out_dir}/shrine", mode: "copy"

    input:
        path main
        each path(summary)

    output:
        path("full_summary.txt", includeInputs: true)

    script:
        """
        cat ${summary} >> ${main}
        echo "Appended ${summary} to ${main}"
        """
}

workflow shrine {
    /*
        Finds optimum DM by strucrure maximization

        Take
            i_ds: path
            	Intensity dynamic spectrum
            tresus: val
            	Time resolution in us
    */
    take:
        i_ds
        tresus

    main:    
        
        // Processing
        generate_profiles(params.label,params.dm_frb,params.dm_low,params.dm_high,params.dm_step,params.dm_count,params.timescale,params.centre_freq_frb,params.bandwidth,i_ds,tresus)
        
        // Select kc (or not)
        if(params.do_uncertainty_min == true){
            minimise_uncertainty(params.label,params.timescale,params.saving,params.force_kc,generate_profiles.out.DMdata,generate_profiles.out.Idata)
            kc = minimise_uncertainty.out.kc
        }
        else{
            kc = params.force_kc
        }

        // Do vary kc. If we didn't minimise uncertainty, run that now because its a useful plot and contains similar information to vary_kc.
        if(params.do_vary_kc == true){
            vary_kc(params.label,params.timescale,params.saving,kc,generate_profiles.out.DMdata,generate_profiles.out.Idata)
            if(params.do_uncertainty_min == false){
                minimise_uncertainty(params.label,params.timescale,params.saving,kc,generate_profiles.out.DMdata,generate_profiles.out.Idata)
            }
        }
        
        smres	= maximise_structure(params.label,params.dm_frb,params.timescale,params.saving,kc,generate_profiles.out.DMdata,generate_profiles.out.Idata)
        
        if(params.do_sn == true){
            maximise_sn(params.label,params.timescale,params.saving,kc,generate_profiles.out.DMdata,generate_profiles.out.Idata)
            if(params.saving == true){
                //If we're doing SN and saving plots, generate this extra plot
                compare_structure_SN(params.label, params.timescale, generate_profiles.out.DMdata, generate_profiles.out.Idata, maximise_sn.out.Ismooth, maximise_structure.out.DM, maximise_sn.out.DM)
            }
        }
        
        qlookplot(i_ds,tresus,smres.resfile)

        //Build the main summary file
        make_summary()

        //Collect summary files from other scripts.
        summaries = generate_profiles.out.summary
        summaries = summaries.concat(maximise_structure.out.summary)
        if(params.do_sn){
            summaries = summaries.concat(maximise_sn.out.summary)
        }
        if(params.do_uncertainty_min || params.do_vary_kc){
            summaries = summaries.concat(minimise_uncertainty.out.summary)
        }
        
        //Append to main summary file
        cat_summaries(make_summary.out, summaries)
}
