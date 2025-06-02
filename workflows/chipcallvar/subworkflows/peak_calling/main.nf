/*
============================================================================
	SUBWORKFLOW: Peak Calling
=============================================================================
*/

include { MACS3_CALLPEAK } from '../../../../modules/local/macs3/callpeak'


workflow PEAK_CALLING {
    take:
    ch_bam // channel: [ meta, bam, bai ]
    
    main:

    // Control samples (input/control)
    ch_control = ch_bam
        .filter { meta, bam, bai -> !meta.control || meta.control == "" }
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }
    
    // Treatment samples (ChIP samples)
    ch_treatment = ch_bam
        .filter { meta, bam, bai -> meta.control && meta.control != "" }
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }

    ch_matched = ch_treatment
        .join(ch_control) // join by patient
        .map { key, treat_meta, treat_bam, treat_bai, ctrl_meta, ctrl_bam, ctrl_bai  ->
            [treat_meta, treat_bam, treat_bai, ctrl_bam, ctrl_bai]
        }    
    // ch_matched.view { "Matched channel: $it" }
    callpeak = MACS3_CALLPEAK(ch_matched)
    
    emit:
    peaks = callpeak.peaks  // [ meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ]
    stats = callpeak.peak_stats
}
