/*
=============================================================================
        MODULE: VCF2MAF
=============================================================================
*/


process VCF2MAF {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/mafs_annotated/${meta.caller}/${meta.patient}", mode: 'copy'
    container "${params.VCF2MAF_CONTAINER}"

    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("${meta.patient}.${meta.caller}.vep.maf"), emit: maf
    
    script:
    """
    # Uncompress the vcf file
    bgzip -c -d ${vcf} > ${meta.patient}.vcf
    
    # Run vcf2maf
    vcf2maf.pl \
        --input-vcf ${meta.patient}.vcf \
        --output-maf ${meta.patient}.${meta.caller}.vep.maf \
        --tumor-id ${meta.patient} \
        --ref-fasta ${params.REFERENCE_GENOME} \
        --vep-data ${params.VEP_CACHE} \
        --vep-path ${params.VEP_PATH} \
        --species homo_sapiens \
        --ncbi-build ${params.ASSEMBLY} \
        --inhibit-vep
    
    # Compress the maf file
    # gzip ${meta.patient}.macs3.vep.maf
    """
}
