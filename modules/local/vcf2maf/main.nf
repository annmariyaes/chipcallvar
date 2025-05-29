// chipcallvar/modules/nf-core/vcf2maf/main.nf

process VCF2MAF {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/mafs_annotated/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.sample}_vep_ann.maf.gz"), emit: maf

    script:
    """
    # Uncompress the vcf file
    bgzip -c -d ${vcf} > ${meta.sample}.vcf
 
    perl /storage/tools/vcf2maf_v1.6.22/mskcc-vcf2maf-f6d0c40/vcf2maf.pl \
        --input-vcf ${meta.sample}.vcf \
        --output-maf ${meta.sample}_vep_ann.maf \
        --tumor-id ${meta.sample} \
        --ref-fasta ${params.REFERENCE_FASTA} \
        --vep-data ${params.VEP_CACHE} \
        --vep-path ${params.VEP_PATH} \
        --ncbi-build GRCh38 \
        --inhibit-vep  # Use this if VEP has already been run

    # Use pigz for parallel gzip compression
    pigz -p ${task.cpus} ${meta.sample}_vep_ann.maf || gzip ${meta.sample}_vep_ann.maf
    """
}
