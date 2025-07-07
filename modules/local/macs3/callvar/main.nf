/*
=============================================================================
        MODULE: MACS3 callvar
=============================================================================
*/


process MACS3_CALLVAR {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/variant_calling/macs3/${meta.id}", mode: 'copy'
    container "${params.MACS3_CONTAINER}"
    
    input:
    tuple val(meta), path(narrow_peaks), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais)
    
    output:
    tuple val(meta), path("${meta.id}.macs3.vcf"), emit: vcf
    
    script:
    def ctrl_flag = ctrl_bams ? "--control ${ctrl_bams.join(' ')}" : ''
    def treat_files = treat_bams.join(' ')
    
    """    
    # Run MACS3 callvar for this chromosome
    macs3 callvar \\
        --peak ${narrow_peaks} \\
        --treatment ${treat_files} \\
        ${ctrl_flag} \\
        --multiple-processing ${task.cpus} \\
        --outdir . \\
        --ofile ${meta.id}.macs3.vcf \\
        --verbose 2 
    """
}
