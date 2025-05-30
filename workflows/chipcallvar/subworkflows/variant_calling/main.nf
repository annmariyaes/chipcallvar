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
    ch_callvar = ch_peaks.filter { meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        peaks.exists() && treat_bams.exists()
    }
    
    // Call variants
    callvar = MACS3_CALLVAR(ch_callvar)
       
    emit:
    vcf = callvar.vcf  // [ meta, annotated_vcf ]
}
