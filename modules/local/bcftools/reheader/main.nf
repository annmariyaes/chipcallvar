/*
=============================================================================
        MODULE: BCFTOOLS 
============================================================================= 
Reheader the VCF file (wrong sample name, missing headers(SB,PL,AD)
Extract SB values and convert them to AD values
Annotate the VCF with these AD values
Calculate VAF values from the AD field
When converting a VCF to a MAF file, the t_ref_count and t_alt_count columns are typically derived from INFO or FORMAT fields that contain read depth (AD, DP, or SB).
*/


process BCFTOOLS_REHEADER {
    tag "$meta.patient"
    publishDir "${params.OUTDIR}/annotation/ensembl_vep/${meta.patient}", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("${meta.patient}.macs3.vep.filled.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.patient}.macs3.vep.filled.vcf.gz.tbi"), emit: vcf_tbi
    
    script:
    """
    # Extract the original VCF header
    bcftools view -h $vcf > header.txt
    
    # Add missing header lines (INFO/SB, FORMAT/AD, FORMAT/PL)
    sed -i '/^##contig=<ID=chr1/i ##INFO=<ID=SB,Number=4,Type=Integer,Description="Strand Bias">' header.txt
    sed -i 's/##FORMAT=<ID=PL,Number=3/##FORMAT=<ID=PL,Number=G/' header.txt
    sed -i '/^##FORMAT=/a ##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Strand bias counts copied from INFO/SB">' header.txt
    sed -i '/^##FORMAT=/a ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depth calculated from INFO/SB">' header.txt
    
    # Create a sample rename mapping file
    SAMPLE=\$(bcftools query -l $vcf)
    sed -i "s/\tSAMPLE\$/\t${meta.patient}/" header.txt
    
    # Reheader the VCF with the updated header and new sample name
    bcftools reheader -h header.txt -o ${meta.patient}.reheader.vcf $vcf
    
    # Move INFO/SB to FORMAT/SB
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SB\\n' ${meta.patient}.reheader.vcf | awk '
    BEGIN { OFS="\\t" }
    {
        if (\$5 != "." && split(\$5, sb, ",") == 4) {
            print \$1, \$2, \$3, \$4, sb[1] "," sb[2] "," sb[3] "," sb[4];
        }
    }' > ${meta.patient}.SB.values.txt
    
    # Compress and index SB file
    singularity exec -B /storage:/storage ${params.HTSLIB_CONTAINER} bgzip ${meta.patient}.SB.values.txt
    singularity exec -B /storage:/storage ${params.HTSLIB_CONTAINER} tabix -s1 -b2 -e2 ${meta.patient}.SB.values.txt.gz
    
    # Extract AD from SB and save to a temporary file
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SB\\n' ${meta.patient}.reheader.vcf | awk '
    BEGIN { OFS="\\t" }
    {
        if (\$5 != "." && split(\$5, sb, ",") == 4) {
            ref_count = sb[1] + sb[2];  # Sum of reference allele reads
            alt_count = sb[3] + sb[4];  # Sum of alternate allele reads
            print \$1, \$2, \$3, \$4, ref_count "," alt_count;
        }
    }' > ${meta.patient}.AD.values.txt
    
    # Compress and index AD file
    singularity exec -B /storage:/storage ${params.HTSLIB_CONTAINER} bgzip ${meta.patient}.AD.values.txt
    singularity exec -B /storage:/storage ${params.HTSLIB_CONTAINER} tabix -s1 -b2 -e2 ${meta.patient}.AD.values.txt.gz
    
    # Annotate VCF with AD, SB fields
    bcftools annotate -a ${meta.patient}.AD.values.txt.gz -c CHROM,POS,REF,ALT,FORMAT/AD -o ${meta.patient}.with.ad.vcf ${meta.patient}.reheader.vcf
    bcftools annotate -a ${meta.patient}.SB.values.txt.gz -c CHROM,POS,REF,ALT,FORMAT/SB -o ${meta.patient}.with.sb.vcf ${meta.patient}.with.ad.vcf
    
    # Create FORMAT/VAF from FORMAT/AD
    bcftools +fill-tags ${meta.patient}.with.sb.vcf -Oz -o ${meta.patient}.macs3.vep.filled.vcf.gz -t FORMAT/VAF
    singularity exec -B /storage:/storage ${params.HTSLIB_CONTAINER} tabix -p vcf ${meta.patient}.macs3.vep.filled.vcf.gz
    
    # Cleanup
    rm header.txt ${meta.patient}.AD.values.txt.gz ${meta.patient}.AD.values.txt.gz.tbi ${meta.patient}.SB.values.txt.gz ${meta.patient}.SB.values.txt.gz.tbi ${meta.patient}.reheader.vcf ${meta.patient}.with.ad.vcf ${meta.patient}.with.sb.vcf
    """
}
