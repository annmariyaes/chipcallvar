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
    tuple val(meta), path(peaks), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais)
    
    output:
    tuple val(meta), path("${meta.id}.macs3.vcf"), emit: vcf
    
    script:
    def ctrl_flag = ctrl_bams ? "--control ${ctrl_bams.join(' ')}" : ''
    def treat_files = treat_bams.join(' ')
    
    """    
    # Run MACS3 callvar for this chromosome
    macs3 callvar \\
        --peak ${peaks} \\
        --treatment ${treat_files} \\
        ${ctrl_flag} \\
        --multiple-processing ${task.cpus} \\
        --outdir . \\
        --ofile ${meta.id}.macs3.vcf \\
        --verbose 2 
    """
}



// Run MACS3 callvar per chromosome
process MACS3_CALLVAR_CHR {
    tag "${meta.id}_${chr}"
    container "${params.MACS3_CONTAINER}"
    label 'process_medium'
    
    input:
    tuple val(meta), val(chr), path(chr_bed), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais)
    
    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    
    script:
    def ctrl_flag = ctrl_bams ? "--control ${ctrl_bams.join(' ')}" : ''
    def treat_files = treat_bams.join(' ')
    def output_name = "${meta.id}_${chr}.macs3.vcf"
    
    """
    # Check if BED file is not empty
    if [[ -s ${chr_bed} ]]; then
        macs3 callvar \\
            --peak ${chr_bed} \\
            --treatment ${treat_files} \\
            ${ctrl_flag} \\
            --multiple-processing ${task.cpus} \\
            --outdir . \\
            --ofile ${output_name} \\
            --verbose 2
    else
        # Create empty VCF if no peaks for this chromosome
        touch ${output_name}
    fi
    """
}
