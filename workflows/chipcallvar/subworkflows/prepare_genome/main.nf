/*
=============================================================================
	SUBWORKFLOW: Reference Preparation
=============================================================================
*/

include { BWAMEM2_INDEX } from '../../../../modules/local/bwamem2/index'


workflow PREPARE_GENOME {
    take:
    ch_reference // channel: path to reference genome
    
    main:
    ch_index = BWAMEM2_INDEX(ch_reference)
    
    emit:
    index = ch_index.index
}