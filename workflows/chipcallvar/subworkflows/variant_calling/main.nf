/*
=============================================================================
    SUBWORKFLOW: Variant Calling
=============================================================================
*/


include { MACS3_CALLVAR } from '../../../../modules/local/macs3/callvar'


workflow VARIANT_CALLING {
    take:
    ch_peaks // channel: [ meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ]
    
    main:
    // Filter valid peaks for variant calling
    ch_callvar = ch_peaks.map { tuple_data ->
            if (tuple_data.size() == 6) {
                // No control sample: [meta, peaks, treat_meta, treat_bam, treat_bai, null]
                def (patient, treat_meta, treat_bam, treat_bai, null_value) = tuple_data
                def updated_meta = treat_meta + [has_control: false]
                [updated_meta, treat_bam, treat_bai, [], []]
            } else {
                // Has control sample: [meta, peaks, treat_meta, treat_bam, treat_bai, ctrl_meta, ctrl_bam, ctrl_bai]
                def (patient, treat_meta, treat_bam, treat_bai, ctrl_meta, ctrl_bam, ctrl_bai) = tuple_data
                def updated_meta = treat_meta + [has_control: true]
                [updated_meta, treat_bam, treat_bai, ctrl_bam, ctrl_bai]
            }
        }.view { "Joined for MACS3 callvar: $it" }
    
    // Call variants
    callvar = MACS3_CALLVAR(ch_callvar)
       
    emit:
    vcf = callvar.vcf  // [ meta, annotated_vcf ]
}
