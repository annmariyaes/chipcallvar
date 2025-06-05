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
    bgzip -c -d ${vcf} > ${meta.patient}.vcf

    perl vcf2maf.pl \
        --input-vcf ${vcf} \
        --output-maf ${meta.patient}.macs3.vep.maf \
        --tumor-id ${meta.patient} \
        --ref-fasta ${params.REFERENCE_GENOME} \
        --vep-data ${params.VEP_CACHE} \
        --vep-path ${params.VEP_PATH} \
        --ncbi-build GRCh38 \
        --inhibit-vep  # Use this if VEP has already been run

    # Use pigz for parallel gzip compression
    pigz -p ${task.cpus} ${meta.patient}.macs3.vep.maf || gzip ${meta.patient}.macs3.vep.maf
    """
}
