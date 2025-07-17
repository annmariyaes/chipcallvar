/*
=============================================================================
        SUBWORKFLOW: Preprocessing
=============================================================================
*/


include { GATK_MARK_DUPLICATES } from '../../../../modules/local/gatk/markduplicates'
include { GATK_BQSR } from '../../../../modules/local/gatk/bqsr'
include { CREATE_INTERVALS_BED } from '../../../../modules/local/bedtools/makewindows'
include { GATK_CREATE_SEQUENCE_DICTONARY } from '../../../../modules/local/gatk/createsequencedictionary'

workflow PRE_PROCESSING {
    take:
    ch_merged
    ch_reference
    ch_fai
    
    main:
    if (!(params.skip_tools && params.skip_tools.split(',').contains('markduplicates') && params.skip_tools.split(',').contains('bqsr'))) {
        GATK_MARK_DUPLICATES(ch_merged)
        ch_dbsnp = Channel.fromPath(params.DBSNP, checkIfExists: true)
        ch_dbindel = Channel.fromPath(params.MILLS_1000G, checkIfExists: true)
        // Combine the channels properly for GATK_BQSR
        ch_bqsr = GATK_MARK_DUPLICATES.out.duplicates_marked
                                      .combine(ch_dbsnp)
                                      .combine(ch_dbindel)
        ch_preprocessed = GATK_BQSR(ch_bqsr)
    }
    else {
        ch_preprocessed = ch_merged
    }
    
    ch_chunks = CREATE_INTERVALS_BED(ch_fai)
    ch_dict = GATK_CREATE_SEQUENCE_DICTONARY(ch_reference)
    
    emit:
    preprocessed = ch_preprocessed
    intervals = ch_chunks.intervals
    dict = ch_dict.dict
}
