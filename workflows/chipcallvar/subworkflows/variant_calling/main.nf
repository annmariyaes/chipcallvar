/*
=============================================================================
    SUBWORKFLOW: Variant Calling
=============================================================================
*/

include { SPLIT_PEAKS_BY_CHR } from '../../../../modules/local/awk'
include { MACS3_CALLVAR_CHR } from '../../../../modules/local/macs3/callvar'
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
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
  
    // Control samples (input/control) - samples that ARE controls
    ch_control = ch_bam
        .filter { meta, bam, bai ->
            meta.containsKey('control') &&
            meta.control != null &&
            meta.control != ""
        }
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }

    // Join treatment and control samples by patient
    ch_joined = ch_treatment
        .join(ch_control, remainder: true)
        .map { tuple_data ->
            if (tuple_data.size() == 5) {
                // No control sample
                def (id, treat_meta, treat_bam, treat_bai, null_value) = tuple_data
                def updated_meta = treat_meta + [has_control: false]
                [updated_meta, treat_bam, treat_bai, [], []]
            }
            else {
                // Has control sample
                def (id, treat_meta, treat_bam, treat_bai, ctrl_meta, ctrl_bam, ctrl_bai) = tuple_data
                def updated_meta = treat_meta + [has_control: true]
                [updated_meta, treat_bam, treat_bai, ctrl_bam, ctrl_bai]
            }
        }

    // Prepare channels for MACS3
    ch_peaks_keyed = ch_peaks.map { meta, peaks ->
        [meta.id, meta, peaks]
    }

    ch_bams_keyed = ch_joined.map { meta, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        [meta.id, meta, treat_bams, treat_bais, ctrl_bams, ctrl_bais]
    }

    ch_callvar = ch_peaks_keyed
        .join(ch_bams_keyed)
        .map { id, peaks_meta, peaks, bams_meta, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
            [peaks_meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais]
        }

    if (params.tools && 'macs3' in params.tools.split(',')) {
        // variant calling without intervals
        if ((params.no_intervals)) {
            MACS3_CALLVAR(ch_callvar)
            ch_macs3_grouped = MACS3_CALLVAR.out.vcf
                .groupTuple(by: 0)
                .map { meta, vcfs -> [meta, vcfs.sort()] }
            
            // Concatenate VCF files if there are multiple
            BCFTOOLS_CONCAT_MACS3(ch_macs3_grouped, 'macs3')
            ch_vcf = ch_vcf.mix(
                BCFTOOLS_CONCAT_MACS3.out.vcf.map { meta, vcf ->
                    [meta + [caller: 'macs3'], vcf]
                })
        }
        // MACS3 with chromosome-level parallelization
        else {
            // Split peaks by chromosome
            SPLIT_PEAKS_BY_CHR(ch_callvar.map { meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
                [meta, peaks]
            })
            // Create channel with chromosome-specific peaks and BAM files
            ch_chr_callvar = SPLIT_PEAKS_BY_CHR.out.chr_peaks
                .transpose()
                .map { meta, chr_peaks_file ->
                    def chr = chr_peaks_file.name.replaceAll(/.*_chr(.+)\.bed/, '$1')
                    [meta.id, meta, chr, chr_peaks_file]
                }
                .combine(ch_bams_keyed.map { id, meta, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
                    [id, treat_bams, treat_bais, ctrl_bams, ctrl_bais]
                }, by: 0)  // Combine by first element (id)
                .map { id, meta, chr, chr_peaks_file, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
                    [meta, chr, chr_peaks_file, treat_bams, treat_bais, ctrl_bams, ctrl_bais]
            }
            // Run MACS3 CALLVAR for each chromosome
            MACS3_CALLVAR_CHR(ch_chr_callvar)
            // Group VCF files by sample for concatenation
            ch_macs3_grouped = MACS3_CALLVAR_CHR.out.vcf
                .groupTuple(by: 0) // Group by meta
                .map { meta, vcfs ->
                    // Sort VCF files by chromosome order for proper concatenation
                    def sorted_vcfs = vcfs.sort { a, b ->
                        def chr_a = a.name.replaceAll(/.*_chr(.+)\.macs3\.vcf/, '$1')
                        def chr_b = b.name.replaceAll(/.*_chr(.+)\.macs3\.vcf/, '$1')
                        // Handle numeric vs non-numeric chromosome names
                        if (chr_a.isNumber() && chr_b.isNumber()) {
                            return chr_a.toInteger() <=> chr_b.toInteger()
                        } else if (chr_a.isNumber() && !chr_b.isNumber()) {
                            return -1
                        } else if (!chr_a.isNumber() && chr_b.isNumber()) {
                            return 1
                        } else {
                            return chr_a <=> chr_b
                        }
                    }
                    [meta, sorted_vcfs]
                }
            
            // Concatenate chromosome-specific VCF files
            BCFTOOLS_CONCAT_MACS3(ch_macs3_grouped, 'macs3')
            ch_vcf = ch_vcf.mix(
                BCFTOOLS_CONCAT_MACS3.out.vcf.map { meta, vcf ->
                    [meta + [caller: 'macs3'], vcf]
                })
        }
    }


    // Prepare channel for GATK/FreeBayes
    ch_vcall = ch_callvar.map { meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        [meta, peaks, treat_bams, treat_bais]
    }
    ch_indexes = ch_index.combine(ch_dict)
    ch_combined = ch_vcall.combine(ch_indexes)

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
