/*
=============================================================================
        MODULE: SAMTOOLS stats, flagstat
=============================================================================
*/

process SAMTOOLS_STATS {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'
    container "${params.SAMTOOLS_CONTAINER}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_bam_stat.txt"), path("${meta.id}_bam_flagstat.txt"), emit: bam_stats

    script:
    """
    # summary of alignment statistics (total reads, mapped reads, properly paired, duplicates, etc.)
    samtools stats $bam > "${meta.id}_bam_stat.txt"

    samtools flagstat $bam > "${meta.id}_bam_flagstat.txt"
    """
}
