/*
=============================================================================
        SUBWORKFLOW: MAF files
=============================================================================
*/

include { VCF2MAF } from '../../../../modules/local/vcf2maf'
include { MAFTOOLS } from '../../../../modules/local/maftools'

workflow MAF_PROCESSING {
    take:
    ch_vcf // channel: [ meta, vcf ]

    main:
    vcf2maf = VCF2MAF(ch_vcf)

    emit:
    maf = vcf2maf.maf  // [ meta, maf ]
}


