/*
=============================================================================
        SUBWORKFLOW: Preprocessing
=============================================================================
*/


include { CREATE_INTERVALS_BED } from '../../../../modules/local/bedtools/makewindows'


workflow PRE_PROCESSING {
    take:
    ch_fai

    main:
    ch_chunks = CREATE_INTERVALS_BED(ch_fai)


    emit:
    intervals = ch_chunks.intervals

}
