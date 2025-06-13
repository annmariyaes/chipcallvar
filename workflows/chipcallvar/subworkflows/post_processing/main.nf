/*
=============================================================================
	SUBWORKFLOW: VCF Post-processing
=============================================================================
*/
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_MACS3 } from '../../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_MUTECT2 } from '../../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_FREEBAYES } from '../../../../modules/local/bcftools/reheader'


workflow VCF_POSTPROCESSING {
    take:
    ch_vep

    main:
    ch_vcf = Channel.empty()

    if (params.tools && params.tools.split(',').contains('macs3')) {
        vep_macs3 = BCFTOOLS_REHEADER_MACS3(ch_vep, 'macs3')
        ch_vcf = ch_vcf.mix(vep_macs3.vcf)
    }

    if (params.tools && params.tools.split(',').contains('mutect2')) {
        vep_mutect2 =BCFTOOLS_REHEADER_MUTECT2(ch_vep, 'mutect2')
        ch_vcf = ch_vcf.mix(vep_mutect2.vcf)
    }

    if (params.tools && params.tools.split(',').contains('freebayes')) {
        vep_freebayes = BCFTOOLS_REHEADER_FREEBAYES(ch_vep, 'freebayes')
        ch_vcf = ch_vcf.mix(vep_freebayes.vcf)
    }

    emit:
    vcf = ch_vcf
}
