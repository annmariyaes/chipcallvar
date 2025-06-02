/*
=============================================================================
        MODULE: BCFTOOLS view, reheader
=============================================================================
*/


// t_ref_count and t_alt_count columns are typically derived from INFO or FORMAT fields that contain read depth (AD, DP, or SB).

process BCFTOOLS {
    tag "$meta.patient"
    publishDir "${params.OUTDIR}/mafs-annotated/${meta.patient}", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.patient}_reheader.vcf"), emit: vcf

    script:
    """
    # Extract the original VCF header
    bcftools view -h $vcf > header.txt

    # Add missing header lines if they are not already present
    grep -q '##INFO=<ID=SB' header.txt || echo '##INFO=<ID=SB,Number=4,Type=Integer,Description="Strand Bias">' >> header.txt
    grep -q '##FORMAT=<ID=PL' header.txt || echo '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">' >> header.txt

    # Create a sample rename mapping file 
    echo "${meta.patient}" > new_sample_name.txt

    # Reheader the VCF with the updated header and new sample name
    bcftools reheader -h header.txt -s new_sample_name.txt $vcf > ${meta.patient}_reheader.vcf

    # Cleanup
    rm header.txt sample_name.txt new_sample_name.txt
    """
}
