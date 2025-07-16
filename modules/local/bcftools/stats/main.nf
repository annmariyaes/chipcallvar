/*
=============================================================================
        MODULE: BCFTOOLS stats
=============================================================================
*/

process BCFTOOLS_STATS {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/multiqc/${meta.caller}", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    val caller

    output:
    tuple val(meta), path("*.vcf.stats.txt"), emit: vcf_stats

    script:
    """
    # Run bcftools stats for MultiQC
    bcftools stats $vcf > "${meta.id}.${meta.replicate}.${meta.caller}.vcf.stats.txt"
    """
}
