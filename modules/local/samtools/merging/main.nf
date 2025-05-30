/*
=============================================================================
        MODULE: SAMTOOLS merge, sort, index
=============================================================================
*/

process SAMTOOLS_MERGE {

    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/preprocessing/merged/${meta.sample}", mode: 'copy'
    container "${params.SAMTOOLS_CONTAINER}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.bam"), path("${meta.sample}.bam.bai"), emit: merged

    script:
    def bams = bam.join(' ')
    """
    echo "Processing BAMs: ${bams}"
    # Merge the technical replicates bam files
    samtools merge -o "${meta.sample}.merged.bam" ${bams}

    # Sort the BAM file
    samtools sort "${meta.sample}.merged.bam" -o "${meta.sample}.merged.sorted.bam"

    # Index the sorted BAM file
    samtools index "${meta.sample}.merged.sorted.bam"

    # Rename sorted BAM as final output
    mv "${meta.sample}.merged.sorted.bam" "${meta.sample}.bam"
    mv "${meta.sample}.merged.sorted.bam.bai" "${meta.sample}.bam.bai"
    """
}
