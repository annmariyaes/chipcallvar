nextflow.enable.dsl=2

include { QUALITY_CONTROL } from '../../workflows/chipcallvar/subworkflows/quality_control'
include { PREPARE_GENOME } from '../../workflows/chipcallvar/subworkflows/prepare_genome'
include { ALIGN_AND_PROCESS } from '../../workflows/chipcallvar/subworkflows/align_process'
include { PRE_PROCESSING } from '../../workflows/chipcallvar/subworkflows/pre_processing'
include { BAM_MERGING } from '../../workflows/chipcallvar/subworkflows/merging'
include { PEAK_CALLING } from '../../workflows/chipcallvar/subworkflows/peak_calling'
include { VARIANT_CALLING } from '../../workflows/chipcallvar/subworkflows/variant_calling'
include { BAM_STATS; VCF_STATS }  from '../../workflows/chipcallvar/subworkflows/stats'
include { VARIANT_ANNOTATION } from '../../workflows/chipcallvar/subworkflows/variant_annotation'
include { VCF_POSTPROCESSING } from '../../workflows/chipcallvar/subworkflows/post_processing'
include { VARIANT_FILTERING } from '../../workflows/chipcallvar/subworkflows/variant_filtering'
include { DOWNSTREAM_ANALYSIS } from '../../workflows/chipcallvar/subworkflows/downstream_analysis'
include { MULTIQC } from '../../modules/local/multiqc'


/*
========================================================================================
   MAIN WORKFLOW - Variant Calling takes fastq as input
========================================================================================
*/

workflow CHIP_SEQ_FASTQ_VARIANT_CALLING {
    take:
    ch_input
    
    main:
    QUALITY_CONTROL(ch_input)
    ch_reference = Channel.fromPath(params.REFERENCE_GENOME, checkIfExists: true)
   
    PREPARE_GENOME(ch_reference)
    ALIGN_AND_PROCESS(ch_input, PREPARE_GENOME.out.index)
    BAM_STATS(ALIGN_AND_PROCESS.out.merged)

    ch_fai = Channel.fromPath(params.GENOME_FAI, checkIfExists: true)
    PRE_PROCESSING(ALIGN_AND_PROCESS.out.merged, ch_reference, ch_fai)

    PEAK_CALLING(PRE_PROCESSING.out.preprocessed)
    VARIANT_CALLING(PRE_PROCESSING.out.preprocessed, PEAK_CALLING.out.peaks, PREPARE_GENOME.out.index, PRE_PROCESSING.out.intervals, PRE_PROCESSING.out.dict)
    VARIANT_ANNOTATION(VARIANT_CALLING.out.vcf)
    VARIANT_FILTERING(VARIANT_ANNOTATION.out.vcf)

    VCF_STATS(VARIANT_FILTERING.out.vcf)

    ch_tpm = Channel.fromPath(params.tpm, checkIfExists: true)
    DOWNSTREAM_ANALYSIS(VARIANT_FILTERING.out.vcf, ch_tpm)
    
    // Combine all channels first
    all_files = QUALITY_CONTROL.out.fastqc_zip
    	.mix(BAM_STATS.out.bam_stats1)
    	.mix(BAM_STATS.out.bam_stats2)
    	.mix(VCF_STATS.out.vcf_stats)
    	.mix(VARIANT_ANNOTATION.out.vep_stats)
    	.map { it[1] }
    	.collect()    
    ch_multiqc_config = Channel.fromPath("${workflow.projectDir}/multiqc_config.yaml", checkIfExists: true)

    MULTIQC(all_files, ch_multiqc_config)

    emit:
    vcf_out = VARIANT_FILTERING.out.vcf
    maf_out = DOWNSTREAM_ANALYSIS.out.maf
    multiqc_html = MULTIQC.out.html
}




/*
========================================================================================
   MAIN WORKFLOW - Variant Calling takes fastq as input
========================================================================================
*/

workflow CHIP_SEQ_BAM_VARIANT_CALLING {
    take:
    ch_input
    
    main:
    QUALITY_CONTROL(ch_input)
    ch_reference = Channel.fromPath(params.REFERENCE_GENOME, checkIfExists: true)

    BAM_MERGING(ch_input)
    BAM_STATS(ALIGN_AND_PROCESS.out.merged)
    PEAK_CALLING(ALIGN_AND_PROCESS.out.merged)
    VARIANT_CALLING(PEAK_CALLING.out.peaks, PREPARE_GENOME.out.index, PRE_PROCESSING.out.intervals)
    VARIANT_ANNOTATION(VARIANT_CALLING.out.vcf)
    VARIANT_FILTERING(VARIANT_ANNOTATION.out.vcf)

    VCF_STATS(VARIANT_FILTERING.out.vcf)
    DOWNSTREAM_ANALYSIS(VARIANT_FILTERING.out.vcf)

    ch_multiqc_config = Channel.fromPath("${workflow.projectDir}/multiqc_config.yaml", checkIfExists: true)
    MULTIQC(
        QUALITY_CONTROL.out.fastqc_zip.map { it[1] }.collect(),
        BAM_STATS.out.bam_stats1.map { it[1] }.collect(),
        BAM_STATS.out.bam_stats2.map { it[1] }.collect(),
        VCF_STATS.out.vcf_stats.map { it[1] }.collect(),
        VARIANT_ANNOTATION.out.vep_stats.map { it[1] }.collect(),
        ch_multiqc_config
        )

    emit:
    vcf_out = VARIANT_FILTERING.out.vcf
    maf_out = DOWNSTREAM_ANALYSIS.out.maf
    multiqc_html = MULTIQC.out.html
}



/*
========================================================================================
   ALTERNATIVE WORKFLOW - Variant Annotation
========================================================================================
*/


workflow CHIP_SEQ_VCF_VARIANT_ANNOTATION {
    take: ch_input

    main:
    VARIANT_ANNOTATION(ch_input)
    VARIANT_FILTERING(VARIANT_ANNOTATION.out.vcf)
    VCF_STATS(VARIANT_FILTERING.out.vcf)
    MAF_PROCESSING(VARIANT_ANNOTATION.out.vcf)

    emit:
    vcf_out = VARIANT_ANNOTATION.out.vcf
    maf_out = DOWNSTREAM_ANALYSIS.out.maf
}
