/*
=============================================================================
        MODULE: GATK MarkDuplicates
// Most tools require duplicates to be tagged in mapped reads to reduce bias.
=============================================================================
*/


process GATK_MARK_DUPLICATES {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/preprocessing/mark_duplicates/${meta.sample}", mode: 'copy'
    container "${params.GATK_CONTAINER}"

    input:
    tuple val(meta), path(treat_bams), path(treat_bais)

    output:
    tuple val(meta), path("${meta.sample}.dedup.bam"), path("${meta.sample}.dedup.bam.bai"), path("${meta.sample}.metrics.txt"), emit: duplicates_marked

    script:
    """
    gatk MarkDuplicates \
        -I "${treat_bams}" \
        -O "${meta.sample}.dedup.bam" \
        -M "${meta.sample}.metrics.txt" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT
    """
}
