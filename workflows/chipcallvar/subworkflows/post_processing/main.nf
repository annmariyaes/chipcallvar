/*
=============================================================================
	SUBWORKFLOW: VCF Post-processing
=============================================================================
*/

include { BCFTOOLS_REHEADER } from '../../../../modules/local/bcftools/reheader'
include { VCF2MAF } from '../../../../modules/local/vcf2maf'


workflow VCF_POSTPROCESSING {
    take:
    ch_vcf // channel: [ meta, vcf ]
    
    main:
    bcftools = BCFTOOLS_REHEADER(ch_vcf)
    
    // vcf2maf = VCF2MAF(bcftools.vcf)
    
    emit:
    vcf = bcftools.vcf // [ meta, processed_vcf ]
    // maf = vcf2maf.maf  // [ meta, maf ]
}
