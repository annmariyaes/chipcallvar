/*
=============================================================================
	SUBWORKFLOW: Read Alignment and Processing
=============================================================================
*/

include { BWAMEM2_MEM } from '../../../../modules/local/bwamem2/mem'
include { SAMTOOLS_MAP } from '../../../../modules/local/samtools/mapping'
include { SAMTOOLS_MERGE } from '../../../../modules/local/samtools/merging'



workflow ALIGN_AND_PROCESS {
    take:
    ch_input     // channel: [ meta, [ reads ] ]
    ch_index     // channel: bwa-mem2 index files

    main:
    // Alignment
    ch_bam = BWAMEM2_MEM(ch_input, ch_index)

    // Process SAM files and get bam & bai files
    ch_map = SAMTOOLS_MAP(ch_bam.sam)

    // Group by sample for technical replicates
    ch_grouped = ch_map.map { meta, bam, bai ->
                    tuple(meta.sample, meta, bam, bai)
        }
        .groupTuple()
        .map { sample, metas, bams, bais ->
                // Create new meta with sample name and combine all files
                def merged_meta = metas[0].clone()
                merged_meta.id = sample
                tuple(merged_meta, bams, bais)
        }
    ch_grouped.view()

    // Merge technical replicates
    ch_merge = SAMTOOLS_MERGE(ch_grouped)

    emit:
    merged = ch_merge.merged  // [ meta, merged_bam, merged_bai ]
}
