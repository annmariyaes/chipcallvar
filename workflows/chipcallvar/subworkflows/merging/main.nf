/*
=============================================================================
	SUBWORKFLOW: Technical replicate merging
=============================================================================
*/

include { SAMTOOLS_MERGE } from '../../../../modules/local/samtools/merging'



workflow BAM_MERGING {
    take:
    ch_input     

    main:

    // Simple grouping by sample
    ch_grouped = ch_input
        .map { meta, bam, bai -> tuple(meta.sample, meta, bam, bai) }
        .groupTuple()
        .map { sample, metas, bams, bais ->
            def merged_meta = metas[0].clone()
            merged_meta.id = sample
            tuple(merged_meta, bams, bais)
        }        
    // Merge technical replicates
    ch_merge = SAMTOOLS_MERGE(ch_grouped)

    emit:
    merged = ch_merge.merged  // [ meta, merged_bam, merged_bai ]
}
