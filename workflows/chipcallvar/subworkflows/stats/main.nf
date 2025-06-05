include { SAMTOOLS_STATS } from '../../../../modules/local/samtools/stats'
include { MOSDEPTH_STATS } from '../../../../modules/local/mosdepth'
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
    MOSDEPTH_STATS(ch_merged)    

    emit:
    bam_stats1 = SAMTOOLS_STATS.out.bam_stats
    bam_stats2 = MOSDEPTH_STATS.out.bam_stats
}


workflow VCF_STATS {
    take:
    ch_vep

    main:
    BCFTOOLS_STATS(ch_vep)
    
    emit:
    vcf_stats = BCFTOOLS_STATS.out.vcf_stats
}
