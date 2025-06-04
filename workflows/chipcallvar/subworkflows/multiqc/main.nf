include { BCFTOOLS_STATS } from '../../../../modules/local/bcftools/stats'


workflow {
    fastqc_out = FASTQC(input_reads)
    samtools_out = SAMTOOLS_STATS(aligned_bams)
    bcftools_out = BCFTOOLS_STATS(filtered_vcfs)
    vep_out = VEP_STATS(annotated_vcfs)

    MULTIQC(
        fastqc_out.collect(),
        samtools_out.collect(),
        bcftools_out.collect(),
        vep_out.collect()
    )
}
