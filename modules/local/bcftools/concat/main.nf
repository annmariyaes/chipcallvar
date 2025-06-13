/*
=============================================================================
        MODULE: BCFTOOLS_CONCAT 
Stitching the vcf files together
=============================================================================
*/


process BCFTOOLS_CONCAT {
    tag "${meta.patient}_${caller}"
    publishDir "${params.OUTDIR}/variant_calling/${caller}/${meta.patient}", mode: 'copy'
    container "${params.BCFTOOLS_CONTAINER}"
    label 'process_low'
    
    input:
    tuple val(meta), path(vcfs)
    val caller
    
    output:
    tuple val(meta), path("${meta.patient}.${caller}.merged.vcf.gz"), emit: vcf
    
    script:
    def vcf_list = vcfs.join('\\n')
    """
    echo "${vcf_list}" > vcf_list.txt
    bcftools concat \
        -f vcf_list.txt \
        -O z \
        -o ${meta.patient}.${caller}.merged.vcf.gz
    
    bcftools index ${meta.patient}.${caller}.merged.vcf.gz
    """
}
