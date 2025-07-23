


process SPLIT_PEAKS_BY_CHR {
    tag "${meta.id}"
    container "${params.MACS3_CONTAINER}"
    label 'process_low'
    
    input:
    tuple val(meta), path(peaks)
    
    output:
    tuple val(meta), path("${meta.id}_chr*.bed"), emit: chr_peaks
    
    script:
    """
    # Split peaks by chromosome and create separate BED files
    awk '{print > "${meta.id}_chr" \$1 ".bed"}' ${peaks}
    
    # Remove empty files (chromosomes with no peaks)
    find . -name "${meta.id}_chr*.bed" -empty -delete
    """
}
