nextflow.enable.dsl=2

// Include nessesary scripts
include {do_ref_correlation; do_correlation; get_start_mjd} from '../pipelines/correlate'


// Cards and FPGAs to be processed. Override these in a config file to cut out
// data. The lowest card-fpga pair is used as a reference correlation.
params.cards = ["1", "2", "3", "4", "5", "6", "7"]
cards = Channel.fromList(params.cards)
params.fpgas = ["0", "1", "2", "3", "4", "5"]
fpgas = Channel.fromList(params.fpgas)
card_fpgas = cards.combine(fpgas)
    .filter{ !(it[0] == params.cards.min() & it[1] == params.fpgas.min()) }
ref_card_fpga = cards.min().combine(fpgas.min())


// defaults for parameters
params.binconfig = ''
params.polyco = ''
params.data = ''
// params.outdir = './output'
params.publishDir = './output'
params.hard_copy = false


process get_inttime {

    input:
        val int_time

    output:
        path "inttime.txt", emit: int_time
    
    script:
        """
        echo $int_time >> inttime.txt
        
        """
}


process save_outputs {
    /*
        save outputs
    */

    input:
        path correlated_data
        path binconfig
        path polyco

    script:
        """
        echo "Saving Outputs"
        
        if [ -d $params.publishDir ]; then
            rm -r $params.publishDir
        fi
        mkdir $params.publishDir

        if [ $params.hard_copy == 'true' ]; then
            cp -r \$(realpath *) $params.publishDir/.
        else
            cp -r * $params.publishDir/.
        fi

        ex="Done"



        """



}


workflow {
    /*

        This nextflow script runs a correlation of the frb data outside of 
        the main CELEBI pipeline and is primarily used to quickly get .difx
        data. 


    */

    // get int_time
    get_inttime(1.3824)


    // Get start mjd
    startmjd = get_start_mjd(params.data_frb)

    // Reference correlation
    ref_correlation = do_ref_correlation(params.label, params.data_frb, params.ra_frb, params.dec_frb, 
                        params.binconfig, params.polyco, get_inttime.out.int_time, startmjd, ref_card_fpga,
                        params.fcm).cx_fy
    
    // Do rest of correlations
    correlated_data = do_correlation(params.label, params.data_frb, params.ra_frb, params.dec_frb, 
                        params.binconfig, params.polyco, get_inttime.out.int_time, startmjd, ref_correlation.combine(card_fpgas),
                        params.fcm).cx_fy

    // Combine correlations
    all_correlations = ref_correlation.concat(correlated_data).collect()

    // save outputs
    save_outputs(all_correlations, params.binconfig, params.polyco)


}
