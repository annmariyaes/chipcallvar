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
    tuple val(meta), path("${meta.id}.bam.stat.txt"), path("${meta.id}.bam.flagstat.txt"), emit: bam_stats

    script:
    """
    # summary of alignment statistics (total reads, mapped reads, properly paired, duplicates, etc.)
    samtools stats $bam > "${meta.id}.bam.stat.txt"

    samtools flagstat $bam > "${meta.id}.bam.flagstat.txt"
    """
}
