/*
=============================================================================
        SUBWORKFLOW: Preprocessing
=============================================================================
*/


include { CREATE_INTERVALS_BED } from '../../../../modules/local/bedtools/makewindows'
include { CREATE_SEQUENCE_DICTONARY } from '../../../../modules/local/gatk/createsequencedictionary'


workflow PRE_PROCESSING {
    take:
    ch_reference
    ch_fai

    main:
    ch_chunks = CREATE_INTERVALS_BED(ch_fai)
    ch_dict = CREATE_SEQUENCE_DICTONARY(ch_reference)

    emit:
    intervals = ch_chunks.intervals
    dict = ch_dict.dict
}
