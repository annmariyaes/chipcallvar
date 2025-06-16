/*
=============================================================================
    SUBWORKFLOW: Variant Calling
=============================================================================
*/


include { MACS3_CALLVAR } from '../../../../modules/local/macs3/callvar'
include { GATK_MUTECT2 } from '../../../../modules/local/gatk/mutect2'
include { FREEBAYES } from '../../../../modules/local/freebayes'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_MACS3 } from '../../../../modules/local/bcftools/concat'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_MUTECT2 } from '../../../../modules/local/bcftools/concat'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_FREEBAYES } from '../../../../modules/local/bcftools/concat'



workflow VARIANT_CALLING {
    take:
    ch_bam
    ch_peaks
    ch_index
    ch_interval
    ch_dict

    main:
    ch_vcf = Channel.empty()

    ch_treatment = ch_bam
        .filter { meta, bam, bai ->
            !meta.containsKey('control') ||
            meta.control == null ||
            meta.control == ""
        }
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }

    // Control samples (input/control) - samples that ARE controls
    ch_control = ch_bam
        .filter { meta, bam, bai ->
            meta.containsKey('control') &&
            meta.control != null &&
            meta.control != ""
        }
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }

    // Join treatment and control samples by patient
    ch_joined = ch_treatment
        .join(ch_control, remainder: true)
        .map { tuple_data ->
            if (tuple_data.size() == 5) {
                // No control sample: [patient, treat_meta, treat_bam, treat_bai, null]
                def (patient, treat_meta, treat_bam, treat_bai, null_value) = tuple_data
                def updated_meta = treat_meta + [has_control: false]
                [updated_meta, treat_bam, treat_bai, [], []]
            }
            else {
                // Has control sample: [patient, treat_meta, treat_bam, treat_bai, ctrl_meta, ctrl_bam, ctrl_bai]
                def (patient, treat_meta, treat_bam, treat_bai, ctrl_meta, ctrl_bam, ctrl_bai) = tuple_data
                def updated_meta = treat_meta + [has_control: true]
                [updated_meta, treat_bam, treat_bai, ctrl_bam, ctrl_bai]
            }
        }

    // Prepare channel for MACS3 as require peaks
    ch_peaks_keyed = ch_peaks.map { meta, peaks ->
        [meta.patient, meta, peaks]
    }

    ch_bams_keyed = ch_joined.map { meta, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        [meta.patient, meta, treat_bams, treat_bais, ctrl_bams, ctrl_bais]
    }

    ch_callvar = ch_peaks_keyed
        .join(ch_bams_keyed)
        .map { patient, peaks_meta, peaks, bams_meta, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
            [peaks_meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais]
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
            })
    }


    // Prepare channel for GATK/FreeBayes
    ch_vcall = ch_callvar.map { meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        [meta, peaks, treat_bams, treat_bais]
    }
    ch_indexes = ch_index.combine(ch_dict)
    ch_combined = ch_vcall.combine(ch_indexes).combine(ch_interval)

    // MUTECT2 - Tag with caller name
    if (params.tools && params.tools.split(',').contains('mutect2')) {
        GATK_MUTECT2(ch_combined)
        ch_mutect2_grouped = GATK_MUTECT2.out.vcf
            .groupTuple(by: 0)
            .map { meta, vcfs -> [meta, vcfs.sort()] }

        ch_vcf = ch_vcf.mix(
            GATK_MUTECT2.out.vcf.map { meta, vcf -> 
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

        ch_vcf = ch_vcf.mix(
            FREEBAYES.out.vcf.map { meta, vcf -> 
                [meta + [caller: 'freebayes'], vcf] 
            }
        )
    }

    emit:
    vcf = ch_vcf
}
