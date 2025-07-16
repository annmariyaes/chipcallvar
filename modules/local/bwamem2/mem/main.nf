/*
=============================================================================
        MODULE: BWA-MEM2 Alignment
=============================================================================
*/


process BWAMEM2_MEM {
    tag "${meta.unique_id}"
    publishDir "${params.OUTDIR}/preprocessing/alignment/${meta.unique_id}", mode: 'copy'    
    container "${params.BWAMEM2_CONTAINER}"

    input:
    tuple val(meta), path(fastqs), path(reference), path("*")

    output:
    tuple val(meta), path("${meta.unique_id}.sam"), emit: sam
    
    script:
    """
    # BWA (Burrows-Wheeler Aligner) to align all of the reads to a genome.
    # Convert fastqs into bash array
    fastqs=( ${fastqs} )
    
    if [[ \${#fastqs[@]} -eq 1 ]]; then
        echo "Single-end detected"
        bwa-mem2 mem -t ${task.cpus} -R "@RG\\tID:${meta.unique_id}\\tSM:${meta.sample}" "${reference}" "\${fastqs[0]}" > "${meta.unique_id}.sam"
    else
        echo "Paired-end detected"
        bwa-mem2 mem -t ${task.cpus} -R "@RG\\tID:${meta.unique_id}\\tSM:${meta.sample}" "${reference}" "\${fastqs[0]}" "\${fastqs[1]}" > "${meta.unique_id}.sam"
    fi
    """
}
