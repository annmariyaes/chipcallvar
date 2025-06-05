/*
=============================================================================
        MODULE: VCF2MAF
=============================================================================
*/


process VCF2MAF {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/mafs_annotated/${meta.patient}", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("${meta.patient}.macs3.vep.maf.gz"), emit: maf
    
    script:
    """
    # Uncompress the vcf file
    singularity exec -B /storage:/storage ${params.HTSLIB_CONTAINER} bgzip -c -d ${vcf} > ${meta.patient}.vcf
    
    # Run vcf2maf
    perl /storage/tools/vcf2maf_v1.6.22/mskcc-vcf2maf-f6d0c40/vcf2maf.pl \
        --input-vcf ${meta.patient}.vcf \
        --output-maf ${meta.patient}.macs3.vep.maf \
        --tumor-id ${meta.patient} \
        --ref-fasta ${params.REFERENCE_GENOME} \
        --vep-data ${params.VEP_CACHE} \
        --vep-path ${params.VEP_PATH} \
        --species homo_sapiens \
        --ncbi-build GRCh38 \
        --samtools-exec ${params.HTSLIB_CONTAINER} \
        --tabix-exec ${params.HTSLIB_CONTAINER} \
        --inhibit-vep
    
    # Compress the maf file
    gzip ${meta.patient}.macs3.vep.maf
    """
}
