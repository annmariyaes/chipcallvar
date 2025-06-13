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

    // Branch the channel by caller using the meta.caller field
    ch_vep.branch {
        macs3: it[0].caller == 'macs3'
        mutect2: it[0].caller == 'mutect2'
        freebayes: it[0].caller == 'freebayes'
        other: true
    }.set { ch_branched }

    if (params.tools && params.tools.split(',').contains('macs3')) {
        vep_macs3 = BCFTOOLS_REHEADER_MACS3(ch_branched.macs3, 'macs3')
        ch_vcf = ch_vcf.mix(vep_macs3.vcf)
    }

    emit:
    vcf = ch_vcf
}
