/*
=============================================================================
	SUBWORKFLOW: VCF Post-processing
=============================================================================
*/

include { BCFTOOLS } from '../../../../modules/local/bcftools'
include { VCF2MAF } from '../../../../modules/local/vcf2maf'


workflow VCF_POSTPROCESSING {
    take:
    ch_vcf // channel: [ meta, vcf ]
    
    main:
    bcftools = BCFTOOLS(ch_vcf)
    
    vcf2maf = VCF2MAF(bcftools.vcf)
    
    emit:
    maf = vcf2maf.maf  // [ meta, maf ]
    vcf = bcftools.vcf // [ meta, processed_vcf ]
}
