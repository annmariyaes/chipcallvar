/*
=============================================================================
        MODULE: BCFTOOLS stats
=============================================================================
*/

process BCFTOOLS_STATS {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/multiqc/${meta.patient}", mode: 'copy'
    container "${params.BCFTOOLS_CONTAINER}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.patient}_vcf_stats.txt"), path("${meta.patient}_stats.vchk"), path("outdir"), emit: vcf_stats

    script:
    """
    # Run bcftools stats for MultiQC
    bcftools stats $vcf > "${meta.patient}_vcf_stats.txt"

    # Generate full variant stats for plotting
    bcftools stats -s - $vcf > "${meta.patient}_stats.vchk"

    # Create visual plots
    plot-vcfstats -p outdir "${meta.patient}_stats.vchk"
    """
}
