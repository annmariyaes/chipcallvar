/*
=============================================================================
        MODULE: MOSDEPTH
=============================================================================
*/


process MOSDEPTH_STATS {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'
    container "${params.MOSDEPTH_CONTAINER}"
    
    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: bam_stats
    path "*.mosdepth.summary.txt"
    path "*.mosdepth.region.dist.txt", optional: true
    
    script:
    def prefix = meta.id
    """
    mosdepth ${prefix} ${bam}
    """
}
