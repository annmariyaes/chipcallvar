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
    tuple val(meta), path("${meta.unique_id}.bam"), path("${meta.unique_id}.bam.bai"), emit: mapped
    
    script:
    """
    # Convert SAM to BAM, sort, and output directly to final filename
    samtools view --threads ${task.cpus} -Sb "${sam}" | samtools sort --threads ${task.cpus} -o "${meta.unique_id}.bam" 

    # Index the final BAM file
    samtools index --threads ${task.cpus} "${meta.unique_id}.bam"
    
    # Clean up the original SAM file
    rm -f "${sam}"
    """
}
