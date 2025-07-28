/*
=============================================================================
        MODULE: BCFTOOLS_CONCAT 
Stitching the vcf files together
=============================================================================
*/


process BCFTOOLS_CONCAT {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/variant_calling/${caller}/${meta.id}", mode: 'copy'
    container "${params.BCFTOOLS_CONTAINER}"
    label 'process_low'
    
    input:
    tuple val(meta), path(vcfs)
    val caller
    
    output:
    tuple val(meta), path("${meta.id}.${caller}.merged.vcf.gz"), emit: vcf
    
    script:

    // Read file names from FILE, one file name per line
    def vcf_list = vcfs.join('\n')
    """
    echo "${vcf_list}" > vcf_list.txt
    bcftools concat \
        -f vcf_list.txt \
        -O z \
        -o ${meta.id}.${caller}.merged.vcf.gz
    
    bcftools index ${meta.id}.${caller}.merged.vcf.gz
    """
}
