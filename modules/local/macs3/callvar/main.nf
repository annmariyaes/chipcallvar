/*
=============================================================================
        MODULE: MACS3 callvar
=============================================================================
*/


process MACS3_CALLVAR {
    tag "${meta.patient}_${chr}"
    publishDir "${params.OUTDIR}/variant_calling/macs3/${meta.patient}", mode: 'copy'
    container "${params.MACS3_CONTAINER}"
    
    input:
    tuple val(meta), val(chr), path(chr_peaks), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais)
    
    output:
    tuple val(meta), val(chr), path("${meta.patient}_${chr}.macs3.vcf"), emit: vcf
    tuple val(meta), val(chr), path("${meta.patient}_${chr}.macs3.log"), emit: log
    
    script:
    def ctrl_flag = ctrl_bams ? "--control ${ctrl_bams.join(' ')}" : ''
    def treat_files = treat_bams.join(' ')
    
    """
    # Validate input files
    if [ ! -s ${chr_peaks} ]; then
        echo "No peaks found for chromosome ${chr}" > ${meta.patient}_${chr}.macs3.log
        touch ${meta.patient}_${chr}.macs3.vcf
        exit 0
    fi
    
    # Run MACS3 callvar for this chromosome
    macs3 callvar \\
        --peak ${chr_peaks} \\
        --treatment ${treat_files} \\
        ${ctrl_flag} \\
        --multiple-processing ${task.cpus} \\
        --outdir . \\
        --ofile ${meta.patient}_${chr}.macs3 \\
        --verbose 2 > ${meta.patient}_${chr}.macs3.log 2>&1
    
    # Ensure VCF file exists (MACS3 might create .vcf extension automatically)
    if [ -f "${meta.patient}_${chr}.macs3" ]; then
        mv "${meta.patient}_${chr}.macs3" "${meta.patient}_${chr}.macs3.vcf"
    elif [ ! -f "${meta.patient}_${chr}.macs3.vcf" ]; then
        touch "${meta.patient}_${chr}.macs3.vcf"
    fi
    """
}
