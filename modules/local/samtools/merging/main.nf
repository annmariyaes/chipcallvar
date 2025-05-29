process SAMTOOLS_MERGE {

    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/preprocessing/merged/${meta.sample}", mode: 'copy'
    container "${params.SAMTOOLS_CONTAINER}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.bam"), path("${meta.sample}.bam.bai"), emit: merged

    script:
    """
    samtools merge -r "${meta.sample}.merged.bam" $bam

    # Sort the BAM file
    samtools sort "${meta.sample}.merged.bam" -o "${meta.sample}.merged.sorted.bam"

    # Index the sorted BAM file
    samtools index "${meta.sample}.merged.sorted.bam"

    # Rename sorted BAM as final output
    mv "${meta.sample}.merged.sorted.bam" "${meta.sample}.bam"
    mv "${meta.sample}.merged.sorted.bam.bai" "${meta.sample}.bam.bai"

    # Optionally remove the original SAM file (if you want to clean up)
    rm -f "${sam}"
    """
}
