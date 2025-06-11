/*
=============================================================================
    SUBWORKFLOW: Variant Calling
=============================================================================
*/


include { MACS3_CALLVAR } from '../../../../modules/local/macs3/callvar'
include { GATK_MUTECT2 } from '../../../../modules/local/mutect2'
include { FREEBAYES } from '../../../../modules/local/freebayes'


workflow VARIANT_CALLING {
    take:
    ch_peaks // channel: [ meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ]
    
    main:

    ch_vcf = Channel.empty()

    // // Prepare channel for MACS3 as require peaks
    ch_callvar = ch_peaks.filter { meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        peaks.exists() && treat_bams.exists()
    }

    // Prepare channel for GATK/FreeBayes (they don't need peaks)
    ch_vcall = ch_peaks.map { meta, peaks, treat_bams, treat_bais, ctrl_bams, ctrl_bais ->
        [meta, treat_bams, treat_bais]
    }.filter { meta, treat_bams, treat_bais ->
        treat_bams.exists()
    }
    
    // MACS3
    if (params.tools && params.tools.split(',').contains('macs3')) {
        MACS3_CALLVAR(ch_callvar)
        ch_vcf = ch_vcf.mix(MACS3_CALLVAR.out.vcf)
    }

    // MUTECT2    
    if (params.tools && params.tools.split(',').contains('mutect2')) {
        GATK_MUTECT2(ch_vcall)
        ch_vcf = ch_vcf.mix(GATK_MUTECT2.out.vcf)
    }
    
    // FREEBAYES
    if (params.tools && params.tools.split(',').contains('freebayes')) {
        FREEBAYES(ch_vcall)
        ch_vcf = ch_vcf.mix(FREEBAYES.out.vcf)
    }
    
    // Alternative: Run all tools if no specific tools specified
    if (!params.tools) {
        MACS3_CALLVAR(ch_callvar)
        GATK_MUTECT2(ch_vcall)
        FREEBAYES(ch_vcall)
        
        ch_vcf = ch_vcf
            .mix(MACS3_CALLVAR.out.vcf)
            .mix(GATK_MUTECT2.out.vcf)
            .mix(FREEBAYES.out.vcf)
    }
    
    emit:
    vcf = ch_vcf       
}
