/*
=============================================================================
        MODULE: BCFTOOLS stats
=============================================================================
*/

process BCFTOOLS_STATS {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/multiqc/${caller}", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    val caller

    output:
    tuple val(meta), path("*.vcf.stats.txt"), emit: vcf_stats

    script:
    """
    # Run bcftools stats for MultiQC
    bcftools stats $vcf > "${meta.patient}.${caller}.vcf.stats.txt"
    """
}
