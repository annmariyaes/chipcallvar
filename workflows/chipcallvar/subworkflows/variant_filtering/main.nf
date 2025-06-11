/*
=============================================================================
        SUBWORKFLOW: Variant filtering
=============================================================================
*/

include { BCFTOOLS_PLUGINS } from '../../../../modules/local/bcftools/filtering'


workflow VCF_POSTPROCESSING {
    take:
    ch_vcf // channel: [ meta, vcf ]

    main:
    bcftools = BCFTOOLS_PLUGINS(ch_vcf)

    emit:
    vcf = bcftools.vcf // [ meta, processed_vcf ]
}
