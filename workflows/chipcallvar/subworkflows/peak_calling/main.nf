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
        .map { meta, bam, bai -> [meta.sample, meta, bam, bai] }
    
    // Treatment samples (ChIP samples)
    ch_treatment = ch_bam
        .filter { meta, bam, bai -> meta.control && meta.control != "" }
        .map { meta, bam, bai ->
            [meta.control, meta, bam, bai]  
        }

    ch_control_by_patient = ch_control.map { meta, bam, bai -> tuple(meta.patient, bam, bai) }
    ch_treatment_by_patient = ch_treatment.map { meta, bam, bai -> tuple(meta.patient, meta, bam, bai) }

    ch_matched = ch_treatment_by_patient.join(ch_control_by_patient)
        .map { patient, meta, treat_bam, treat_bai, ctrl_bam, ctrl_bai ->
            tuple(meta, treat_bam, treat_bai, ctrl_bam, ctrl_bai)
        }
    callpeak = MACS3_CALLPEAK(ch_bam)
    
    emit:
    peaks = callpeak.peaks  // [ meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ]
}