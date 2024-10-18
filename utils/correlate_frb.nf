nextflow.enable.dsl=2

// Include nessesary scripts
include {do_ref_correlation; do_correlation; get_start_mjd} from '../pipelines/correlate'
include {correlate_frb} from '../pipelines/mfimage'


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
params.fcm = ''
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

    // correlate frb

    correlate_frb(params.binconfig, params.polyco, params.fcm)

    // save outputs
    save_outputs(all_correlations, params.binconfig, params.polyco)


}
