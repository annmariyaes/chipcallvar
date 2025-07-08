/*
=============================================================================
        MODULE: SAMTOOLS stats, flagstat
=============================================================================
*/

process SAMTOOLS_STATS {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'
    container "${params.SAMTOOLS_CONTAINER}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.bam.stat.txt"), path("${meta.sample}.bam.flagstat.txt"), emit: bam_stats

    script:
    """
    # summary of alignment statistics (total reads, mapped reads, properly paired, duplicates, etc.)
    samtools stats $bam > "${meta.sample}.bam.stat.txt"

    samtools flagstat $bam > "${meta.sample}.bam.flagstat.txt"
    """
}
