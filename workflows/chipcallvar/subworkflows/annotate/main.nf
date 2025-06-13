/*
=============================================================================
    SUBWORKFLOW: Variant Annotation
=============================================================================
*/


include { ENSEMBL_VEP as ENSEMBL_VEP_MACS3 } from '../../../../modules/local/ensembl-vep'
include { ENSEMBL_VEP as ENSEMBL_VEP_MUTECT2 } from '../../../../modules/local/ensembl-vep'
include { ENSEMBL_VEP as ENSEMBL_VEP_FREEBAYES } from '../../../../modules/local/ensembl-vep'

workflow VARIANT_ANNOTATION {
    take:
    ch_vep // channel: [ meta, vcf ]

    main:
    ch_vcf = Channel.empty()
    ch_vep_stats = Channel.empty()

    if (params.tools && params.tools.split(',').contains('macs3')) {
        vep_macs3 = ENSEMBL_VEP_MACS3(ch_vep, 'macs3')
        ch_vcf = ch_vcf.mix(vep_macs3.vcf)
        ch_vep_stats = ch_vep_stats.mix(vep_macs3.vep_stats)
    }

    if (params.tools && params.tools.split(',').contains('mutect2')) {
        vep_mutect2 = ENSEMBL_VEP_MUTECT2(ch_vep, 'mutect2')
        ch_vcf = ch_vcf.mix(vep_mutect2.vcf)
        ch_vep_stats = ch_vep_stats.mix(vep_mutect2.vep_stats)
    }

    if (params.tools && params.tools.split(',').contains('freebayes')) {
        vep_freebayes = ENSEMBL_VEP_FREEBAYES(ch_vep, 'freebayes')
        ch_vcf = ch_vcf.mix(vep_freebayes.vcf)
        ch_vep_stats = ch_vep_stats.mix(vep_freebayes.vep_stats)
    }

    emit:
    vcf = ch_vcf         // [ meta, annotated_vcf ]
    vep_stats = ch_vep_stats
}
