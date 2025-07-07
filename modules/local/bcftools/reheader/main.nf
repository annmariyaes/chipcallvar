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
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/annotation/${caller}/${meta.id}", mode: 'copy'
    container "${params.BCFTOOLS_CONTAINER}"
    
    input:
    tuple val(meta), path(vcf)
    val caller
    
    output:
    tuple val(meta), path("${meta.id}.${caller}.vep.filled.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.${caller}.vep.filled.vcf.gz.tbi"), emit: vcf_tbi
     
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
    sed -i "s/\tSAMPLE\$/\t${meta.id}/" header.txt
    
    # Reheader the VCF with the updated header and new sample name
    bcftools reheader -h header.txt -o ${meta.id}.reheader.vcf $vcf


    # Move INFO/SB to FORMAT/SB
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SB\\n' ${meta.id}.reheader.vcf | awk '
    BEGIN { OFS="\\t" }
    {
        if (\$5 != "." && split(\$5, sb, ",") == 4) {
            print \$1, \$2, \$3, \$4, sb[1] "," sb[2] "," sb[3] "," sb[4];
        }
    }' > ${meta.id}.SB.values.txt

    # Compress and index SB file
    bgzip ${meta.id}.SB.values.txt
    tabix -s1 -b2 -e2 ${meta.id}.SB.values.txt.gz

    # Extract AD from SB and save to a temporary file
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SB\\n' ${meta.id}.reheader.vcf | awk '
    BEGIN { OFS="\\t" }
    {
        if (\$5 != "." && split(\$5, sb, ",") == 4) {
            ref_count = sb[1] + sb[2];  # Sum of reference allele reads
            alt_count = sb[3] + sb[4];  # Sum of alternate allele reads
            print \$1, \$2, \$3, \$4, ref_count "," alt_count;
        }
    }' > ${meta.id}.AD.values.txt

    # Compress and index AD file
    bgzip ${meta.id}.AD.values.txt
    tabix -s1 -b2 -e2 ${meta.id}.AD.values.txt.gz

    # Annotate VCF with AD, SB fields
    bcftools annotate -a ${meta.id}.AD.values.txt.gz -c CHROM,POS,REF,ALT,FORMAT/AD -o ${meta.id}.with.ad.vcf ${meta.id}.reheader.vcf
    bcftools annotate -a ${meta.id}.SB.values.txt.gz -c CHROM,POS,REF,ALT,FORMAT/SB -o ${meta.id}.with.sb.vcf ${meta.id}.with.ad.vcf

    # Create FORMAT/VAF from FORMAT/AD
    bcftools +fill-tags ${meta.id}.with.sb.vcf -Oz -o ${meta.id}.${caller}.vep.filled.vcf.gz -t FORMAT/VAF
    tabix -p vcf ${meta.id}.${caller}.vep.filled.vcf.gz

    # Cleanup
    rm header.txt ${meta.id}.AD.values.txt.gz ${meta.id}.AD.values.txt.gz.tbi ${meta.id}.SB.values.txt.gz ${meta.id}.SB.values.txt.gz.tbi ${meta.id}.reheader.vcf ${meta.id}.with.ad.vcf ${meta.id}.with.sb.vcf    
    """
}

