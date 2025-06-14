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
    
    // Branch the channel by caller using the meta.caller field
    ch_vep.branch {
        macs3: it[0].caller == 'macs3'
        mutect2: it[0].caller == 'mutect2'
        freebayes: it[0].caller == 'freebayes'
        other: true
    }.set { ch_branched }

    if (params.tools && params.tools.split(',').contains('macs3')) {
        stats_macs3 = BCFTOOLS_STATS_MACS3(ch_branched.macs3, 'macs3')
        ch_vcf_stats = ch_vcf_stats.mix(stats_macs3.vcf_stats)
    }

    if (params.tools && params.tools.split(',').contains('mutect2')) {
        stats_mutect2 = BCFTOOLS_STATS_MUTECT2(ch_branched.mutect2, 'mutect2')
        ch_vcf_stats = ch_vcf_stats.mix(stats_mutect2.vcf_stats)
    }

    if (params.tools && params.tools.split(',').contains('freebayes')) {
        stats_freebayes = BCFTOOLS_STATS_FREEBAYES(ch_branched.freebayes, 'freebayes')
        ch_vcf_stats = ch_vcf_stats.mix(stats_freebayes.vcf_stats)
    }

    emit:
    vcf_stats = ch_vcf_stats
}
