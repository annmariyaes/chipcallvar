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
    // Treatment samples (ChIP samples) - samples without control specified
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
            } else {
                // Has control sample: [patient, treat_meta, treat_bam, treat_bai, ctrl_meta, ctrl_bam, ctrl_bai]
                def (patient, treat_meta, treat_bam, treat_bai, ctrl_meta, ctrl_bam, ctrl_bai) = tuple_data
                def updated_meta = treat_meta + [has_control: true]
                [updated_meta, treat_bam, treat_bai, ctrl_bam, ctrl_bai]
            }
        }.view { "Joined for MACS3 callpeak: $it" }
    
    // Run MACS3 callpeak
    callpeak = MACS3_CALLPEAK(ch_joined)
    
    emit:
    peaks = callpeak.peaks
    stats = callpeak.peak_stats
}
