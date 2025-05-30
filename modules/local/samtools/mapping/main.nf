/*
=============================================================================
        MODULE: SAMTOOLS view, sort, index
=============================================================================
*/

process SAMTOOLS_MAP {

    tag "${meta.id}"
    publishDir "${params.OUTDIR}/preprocessing/mapped/${meta.id}", mode: 'copy'
    container "${params.SAMTOOLS_CONTAINER}"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: mapped
    
    script:
    """
    # Convert SAM to BAM
    samtools view -Sb "${sam}" > "${meta.id}.bam"

    # Sort the BAM file
    samtools sort "${meta.id}.bam" -o "${meta.id}.sorted.bam"

    # Index the sorted BAM file
    samtools index "${meta.id}.sorted.bam"

    # Rename sorted BAM as final output
    mv "${meta.id}.sorted.bam" "${meta.id}.bam"
    mv "${meta.id}.sorted.bam.bai" "${meta.id}.bam.bai"

    # Optionally remove the original SAM file (if you want to clean up)
    rm -f "${sam}"
    """
}
