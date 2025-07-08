/*
=============================================================================
        MODULE: MOSDEPTH
=============================================================================
*/


process MOSDEPTH_STATS {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'
    container "${params.MOSDEPTH_CONTAINER}"
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: bam_stats
    path "*.mosdepth.summary.txt"
    path "*.mosdepth.region.dist.txt", optional: true
    
    script:
    // <prefix>: base name for output files (no extension)
    def prefix = meta.sample
    """
    mosdepth ${prefix} ${bam}
    """
}
