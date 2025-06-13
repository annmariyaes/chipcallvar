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

    r_script = file("${projectDir}/bin/variant_calling.R")
    // MAFTOOLS(ch_maf_dir, r_script)
    
    emit:
    maf = vcf2maf.maf  // [ meta, maf ]
}


