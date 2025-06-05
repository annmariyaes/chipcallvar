/*
=============================================================================
        MODULE: VCF2MAF
=============================================================================
*/


process VCF2MAF {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/mafs_annotated/${meta.patient}", mode: 'copy'
    conda 'bioconda::htslib bioconda::perl bioconda::pigz bioconda::samtools'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("${meta.patient}.macs3.vep.maf.gz"), emit: maf
    
    script:
    """
    # Uncompress the vcf file
    gunzip -c ${vcf} > ${meta.patient}.vcf
    
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
        --inhibit-vep
    
    # Compress the output
    pigz -p ${task.cpus} ${meta.patient}.macs3.vep.maf || gzip ${meta.patient}.macs3.vep.maf
    """
}
