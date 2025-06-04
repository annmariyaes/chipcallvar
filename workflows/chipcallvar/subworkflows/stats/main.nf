include { SAMTOOLS_STATS } from '../../../../modules/local/samtools/stats'
include { BCFTOOLS_STATS } from '../../../../modules/local/bcftools/stats'


/*
=============================================================================
        SUBWORKFLOW: BAM & VCF statistics
=============================================================================
*/



workflow BAM_STATS {
    take:
    ch_merged

    main:
    SAMTOOLS_STATS(ch_merged)

}


workflow VCF_STATS {
    take:
    ch_vep

    main:
    BCFTOOLS_STATS(ch_vep)

}
