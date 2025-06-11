/*
=============================================================================
        SUBWORKFLOW: BCFTOOLS variant filtering
=============================================================================
*/

include { BCFTOOLS_PLUGINS } from '../../../../modules/local/bcftools/filtering'


workflow VARIANT_FILTERING {
    take:
    ch_vcf // channel: [ meta, vcf ]

    main:
    bcftools = BCFTOOLS_PLUGINS(ch_vcf)

    emit:
    vcf = bcftools.vcf // [ meta, processed_vcf ]
}
