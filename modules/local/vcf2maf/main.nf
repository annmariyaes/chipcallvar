/*
=============================================================================
        MODULE: VCF2MAF
=============================================================================
*/


process VCF2MAF {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/mafs_annotated/${meta.caller}/${meta.id}", mode: 'copy'
    container "${params.VCF2MAF_CONTAINER}"

    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.vep.maf"), emit: maf
    
    script:
    """
    # Uncompress the vcf file
    bgzip -dc ${vcf} > ${meta.id}.vcf
    
    # Run vcf2maf
    vcf2maf.pl \
        --input-vcf ${meta.id}.vcf \
        --output-maf ${meta.id}.${meta.caller}.vep.maf \
        --tumor-id ${meta.id} \
        --ref-fasta ${params.REFERENCE_GENOME} \
        --vep-data ${params.VEP_CACHE} \
        --cache-version ${params.VEP_VERSION} \
        --vep-path ${params.VEP_PATH} \
        --species homo_sapiens \
        --ncbi-build ${params.ASSEMBLY} \
        --inhibit-vep
    
    # Compress the maf file
    gzip ${meta.id}.${meta.caller}.vep.maf
    """
}
