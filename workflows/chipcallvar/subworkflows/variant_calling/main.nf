/*
=============================================================================
    SUBWORKFLOW: Variant Calling
=============================================================================
*/


include { MACS3_CALLVAR } from '../../../../modules/local/macs3/callvar'
include { GATK_MUTECT2 } from '../../../../modules/local/mutect2'
include { FREEBAYES } from '../../../../modules/local/freebayes'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_MACS3 } from '../../../../modules/local/bcftools/concat'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_MUTECT2 } from '../../../../modules/local/bcftools/concat'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_FREEBAYES } from '../../../../modules/local/bcftools/concat'



workflow VARIANT_CALLING {
    take:
    ch_peaks
    ch_index
    ch_interval

    main:
    ch_vcf = Channel.empty()

    // Prepare channel for MACS3 as require peaks
    ch_callvar = ch_peaks.filter { meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        peaks.exists() && treat_bams.exists()
    }
    ch_combinedd = ch_callvar.combine(ch_interval)

    // MACS3 - Tag with caller name
    if (params.tools && params.tools.split(',').contains('macs3')) {
        MACS3_CALLVAR(ch_combinedd)
        ch_macs3_grouped = MACS3_CALLVAR.out.vcf
            .groupTuple(by: 0)
            .map { meta, vcfs -> [meta, vcfs.sort()] }

        BCFTOOLS_CONCAT_MACS3(ch_macs3_grouped, 'macs3')
        ch_vcf = ch_vcf.mix(
            BCFTOOLS_CONCAT_MACS3.out.vcf.map { meta, vcf -> 
                [meta + [caller: 'macs3'], vcf] 
            }
        )
    }

    // Prepare channel for GATK/FreeBayes
    ch_vcall = ch_peaks.map { meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        [meta, treat_bams, treat_bais]
    }
    ch_dict = Channel.fromPath(params.GENOME_DICT, checkIfExists: true)
    ch_indexes = ch_index.combine(ch_dict)
    ch_combined = ch_vcall.combine(ch_indexes).combine(ch_interval)

    // MUTECT2 - Tag with caller name
    if (params.tools && params.tools.split(',').contains('mutect2')) {
        GATK_MUTECT2(ch_combined)
        ch_mutect2_grouped = GATK_MUTECT2.out.vcf
            .groupTuple(by: 0)
            .map { meta, vcfs -> [meta, vcfs.sort()] }

        BCFTOOLS_CONCAT_MUTECT2(ch_mutect2_grouped, 'mutect2')
        ch_vcf = ch_vcf.mix(
            BCFTOOLS_CONCAT_MUTECT2.out.vcf.map { meta, vcf -> 
                [meta + [caller: 'mutect2'], vcf] 
            }
        )
    }

    // FREEBAYES - Tag with caller name
    if (params.tools && params.tools.split(',').contains('freebayes')) {
        FREEBAYES(ch_combined)
        ch_freebayes_grouped = FREEBAYES.out.vcf
            .groupTuple(by: 0)
            .map { meta, vcfs -> [meta, vcfs.sort()] }

        BCFTOOLS_CONCAT_FREEBAYES(ch_freebayes_grouped, 'freebayes')
        ch_vcf = ch_vcf.mix(
            BCFTOOLS_CONCAT_FREEBAYES.out.vcf.map { meta, vcf -> 
                [meta + [caller: 'freebayes'], vcf] 
            }
        )
    }

    emit:
    vcf = ch_vcf
}
