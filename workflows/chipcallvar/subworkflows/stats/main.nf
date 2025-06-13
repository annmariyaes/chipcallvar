/*
=============================================================================
        SUBWORKFLOW: BAM & VCF statistics
=============================================================================
*/


include { SAMTOOLS_STATS } from '../../../../modules/local/samtools/stats'
include { MOSDEPTH_STATS } from '../../../../modules/local/mosdepth'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_MACS3 } from '../../../../modules/local/bcftools/stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_MUTECT2 } from '../../../../modules/local/bcftools/stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_FREEBAYES } from '../../../../modules/local/bcftools/stats'


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
    ch_vcf_stats = Channel.empty()

    if (params.tools && params.tools.split(',').contains('macs3')) {
        BCFTOOLS_STATS_MACS3(
            ch_vep.filter { meta, vcf -> meta.variantcaller == 'macs3' || meta.tool == 'macs3' },
            'macs3'
        )
        ch_vcf_stats = ch_vcf_stats.mix(BCFTOOLS_STATS_MACS3.out.vcf_stats)
    }
    
    if (params.tools && params.tools.split(',').contains('mutect2')) {
        BCFTOOLS_STATS_MUTECT2(
            ch_vep.filter { meta, vcf -> meta.variantcaller == 'mutect2' || meta.tool == 'mutect2' },
            'mutect2'
        )
        ch_vcf_stats = ch_vcf_stats.mix(BCFTOOLS_STATS_MUTECT2.out.vcf_stats)
    }
    
    if (params.tools && params.tools.split(',').contains('freebayes')) {
        BCFTOOLS_STATS_FREEBAYES(
            ch_vep.filter { meta, vcf -> meta.variantcaller == 'freebayes' || meta.tool == 'freebayes' },
            'freebayes'
        )
        ch_vcf_stats = ch_vcf_stats.mix(BCFTOOLS_STATS_FREEBAYES.out.vcf_stats)
    }

    emit:
    vcf_stats = ch_vcf_stats
}
