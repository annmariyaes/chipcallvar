/*
=============================================================================
    SUBWORKFLOW: Variant Annotation
=============================================================================
*/


include { ENSEMBL_VEP as ENSEMBL_VEP_MACS3 } from '../../../../modules/local/ensembl-vep'
include { BCFTOOLS_REHEADER } from '../../../../modules/local/bcftools/reheader'
include { ENSEMBL_VEP as ENSEMBL_VEP_MUTECT2 } from '../../../../modules/local/ensembl-vep'
include { ENSEMBL_VEP as ENSEMBL_VEP_FREEBAYES } from '../../../../modules/local/ensembl-vep'


workflow VARIANT_ANNOTATION {
    take:
    ch_vep 

    main:
    ch_vcf = Channel.empty()
    ch_vep_stats = Channel.empty()

    // Branch the channel by caller using the meta.caller field
    ch_vep.branch {
        macs3: it[0].caller == 'macs3'
        mutect2: it[0].caller == 'mutect2'
        freebayes: it[0].caller == 'freebayes'
        other: true // catch any untagged VCFs
    }.set { ch_branched }

    if (params.tools && params.tools.split(',').contains('macs3')) {
        vep_macs3 = ENSEMBL_VEP_MACS3(ch_branched.macs3, 'macs3')
        head = BCFTOOLS_REHEADER(vep_macs3.vcf, 'macs3')
        ch_vcf = ch_vcf.mix(head.vcf)
        ch_vep_stats = ch_vep_stats.mix(vep_macs3.vep_stats)
    }

    if (params.tools && params.tools.split(',').contains('mutect2')) {
        vep_mutect2 = ENSEMBL_VEP_MUTECT2(ch_branched.mutect2, 'mutect2')
        ch_vcf = ch_vcf.mix(vep_mutect2.vcf)
        ch_vep_stats = ch_vep_stats.mix(vep_mutect2.vep_stats)
    }

    if (params.tools && params.tools.split(',').contains('freebayes')) {
        vep_freebayes = ENSEMBL_VEP_FREEBAYES(ch_branched.freebayes, 'freebayes')
        ch_vcf = ch_vcf.mix(vep_freebayes.vcf)
        ch_vep_stats = ch_vep_stats.mix(vep_freebayes.vep_stats)
    }

    emit:
    vcf = ch_vcf         
    vep_stats = ch_vep_stats
}
