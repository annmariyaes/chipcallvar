/*
=============================================================================
        MODULE: BCFTOOLS stats
=============================================================================
*/

process BCFTOOLS_STATS {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.patient}_vcf_stats.txt"), path("outdir"), emit: vcf_stats

    script:
    """
    # Run bcftools stats for MultiQC
    bcftools stats $vcf > "${meta.patient}_vcf_stats.txt"
    """
}
