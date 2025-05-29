// chipcallvar/modules/nf-core/bwa/mem/main.nf

nextflow.enable.dsl=2

process BWAMEM2_MEM {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/preprocessing/alignment/${meta.id}", mode: 'copy'    
    container "${params.BWAMEM2_CONTAINER}"
    label 'process_high'

    input:
    tuple val(meta), path(fastqs)
    tuple path(reference), path(amb), path(ann), path(bwt), path(pac), path(sa)
    
    output:
    tuple val(meta), path("${meta.sample}.sam"), emit: sam
    
    script:
    """
    # BWA (Burrows-Wheeler Aligner) to align all of the reads to a genome.
    # Convert fastqs into bash array
    fastqs=( ${fastqs} )
    
    if [[ \${#fastqs[@]} -eq 1 ]]; then
        echo "Single-end detected"
        bwa-mem2 mem -t 4 -R "@RG\\tID:${meta.id}\\tSM:${meta.sample}" "${reference}" "\${fastqs[0]}" > "${meta.sample}.sam"
    else
        echo "Paired-end detected"
        bwa-mem2 mem -t 4 -R "@RG\\tID:${meta.id}\\tSM:${meta.sample}" "${reference}" "\${fastqs[0]}" "\${fastqs[1]}" > "${meta.sample}.sam"
    fi
    """
}
