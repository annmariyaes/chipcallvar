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
    
    ch_maf_dir = vcf2maf.maf.map { meta, maf_file ->
        [meta, maf_file.parent]
    }
    MAFTOOLS(ch_maf_dir)
    
    emit:
    maf = vcf2maf.maf  // [ meta, maf ]
}


