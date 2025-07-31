/*
=============================================================================
        SUBWORKFLOW: Consensus variant calling 
Refers to the process of identifying high-confidence variants by integrating results from multiple variant callers or replicates. 
This approach helps reduce false positives and improves reliability.
=============================================================================
*/

include { BCFTOOLS_MERGING } from '../../../../modules/local/bcftools/merging'


workflow CONSENSUS_CALLING {
    take:
    ch_vcf     // Channel: [meta, vcf_file]
    ch_vcf_tbi // Channel: [meta, tbi_file]
    
    main:
    if (params.merge_vcfs) {
        
        ch_vcf_grouped = ch_vcf
            .map { meta, vcf ->
                // Assume TBI exists alongside VCF
                def tbi = file("${vcf}.tbi")
                def sample = meta.sample
                [sample, [meta: meta, vcf: vcf, tbi: tbi]]
            }
            .groupTuple(by: 0)
            .map { sample, items ->
                // Remove duplicates and sort by caller for consistency
                def unique_items = items.unique { it.vcf.toString() }
                    .sort { it.meta.id } // Sort by caller name
                
                def vcfs = unique_items.collect { it.vcf }
                def tbis = unique_items.collect { it.tbi }
                def callers = unique_items.collect { it.meta.id }
                
                def meta = [
                    sample: sample,
                    id: "merged_${callers.join('_')}",
                    callers: callers
                ]
                [meta, vcfs, tbis]
            }
            // .view()        
        ch_merged_vcfs = BCFTOOLS_MERGING(ch_vcf_grouped)
    }
    else {
        ch_merged_vcfs = Channel.empty()
    }
    
    emit:
    merged = ch_merged_vcfs.vcf
}
