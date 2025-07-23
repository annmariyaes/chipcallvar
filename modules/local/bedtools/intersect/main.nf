



process SPLIT_PEAKS_BY_INTERVAL {
    tag "${meta.id}_${interval_name}"
    publishDir "${params.OUTDIR}/preprocessing/intervals", mode: 'copy'
    container "${params.BEDTOOLS_CONTAINER}"
    
    input:
    tuple val(meta), path(peaks), path(intervals)
    
    output:
    tuple val(meta), path("*_filtered.narrowPeak"), path(intervals)
    
    script:
    interval_name = intervals.baseName
    """
    bedtools intersect -a ${peaks} -b ${intervals} > ${meta.id}_${interval_name}_filtered.narrowPeak
    """
}
