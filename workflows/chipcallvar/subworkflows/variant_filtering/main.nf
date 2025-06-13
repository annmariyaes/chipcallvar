/*
=============================================================================
        SUBWORKFLOW: Variant filtering
=============================================================================
*/

include { BCFTOOLS as BCFTOOLS_MACS3 } from '../../../../modules/local/bcftools/filtering'
include { BCFTOOLS as BCFTOOLS_MUTECT2 } from '../../../../modules/local/bcftools/filtering'
include { BCFTOOLS as BCFTOOLS_FREEBAYES } from '../../../../modules/local/bcftools/filtering'

workflow VARIANT_FILTERING {
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
        vep_macs3 = BCFTOOLS_MACS3(ch_branched.macs3, 'macs3')
        ch_vcf = ch_vcf.mix(vep_macs3.vcf)
    }

    if (params.tools && params.tools.split(',').contains('mutect2')) {
        vep_mutect2 = BCFTOOLS_MUTECT2(ch_branched.mutect2, 'mutect2')
        ch_vcf = ch_vcf.mix(vep_mutect2.vcf)
    }

    if (params.tools && params.tools.split(',').contains('freebayes')) {
        vep_freebayes = BCFTOOLS_FREEBAYES(ch_branched.freebayes, 'freebayes')
        ch_vcf = ch_vcf.mix(vep_freebayes.vcf)
    }

    emit:
    vcf = ch_vcf         
}
