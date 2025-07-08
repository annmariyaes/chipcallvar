/*
=============================================================================
        MODULE: SAMTOOLS merge, sort, index
=============================================================================
*/


process SAMTOOLS_MERGE {
    tag "${meta.unique_id}" 
    publishDir "${params.OUTDIR}/preprocessing/merged/${meta.sample}", mode: 'copy'
    container "${params.SAMTOOLS_CONTAINER}"
    
    input:
    tuple val(meta), path(bams), path(bais)
    
    output:
    tuple val(meta), path("${meta.sample}.bam"), path("${meta.sample}.bam.bai"), emit: merged
    
    script:
    def bam_count = bams instanceof List ? bams.size() : 1
    def bam_files = bams instanceof List ? bams.join(' ') : bams.toString()

    """
    if [ ${bam_count} -eq 1 ]; then
        # Single BAM: direct sort (skip merge step)
        # echo "Single BAM files: ${bams}"
        samtools sort --threads ${task.cpus} -o "${meta.sample}.bam" $bams
    else
        # Multiple BAMs: merge then sort
        # echo "Multiple BAM files: ${bams}"
        samtools merge --threads ${task.cpus} -f -o "${meta.sample}.merged.bam" $bams
        samtools sort --threads ${task.cpus} -o "${meta.sample}.bam" "${meta.sample}.merged.bam"
        rm "${meta.sample}.merged.bam"
    fi

    samtools index --threads ${task.cpus} "${meta.sample}.bam"
    """
}
