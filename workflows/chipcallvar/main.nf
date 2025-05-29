include { QUALITY_CONTROL } from '../../workflows/chipcallvar/subworkflows/quality_control'
include { PREPARE_GENOME } from '../../workflows/chipcallvar/subworkflows/prepare_genome'
include { ALIGN_AND_PROCESS } from '../../workflows/chipcallvar/subworkflows/align_process'
include { PEAK_CALLING } from '../../workflows/chipcallvar/subworkflows/peak_calling'
include { VARIANT_CALLING } from '../../workflows/chipcallvar/subworkflows/variant_calling'
include { VARIANT_ANNOTATION } from '../../workflows/chipcallvar/subworkflows/annotate'
include { VCF_POSTPROCESSING } from '../../workflows/chipcallvar/subworkflows/post_processing'

/*
========================================================================================
   MAIN WORKFLOW - Variant Calling
========================================================================================
*/

workflow CHIP_SEQ_VARIANT_CALLING {
    take:
    ch_input

    main:
    QUALITY_CONTROL(ch_input)

    ch_reference = Channel.fromPath(params.REFERENCE_GENOME, checkIfExists: true)     // Prepare reference
    PREPARE_GENOME(ch_reference)

    ALIGN_AND_PROCESS(ch_input, PREPARE_GENOME.out.index)

    PEAK_CALLING(ALIGN_AND_PROCESS.out.merged)

    VARIANT_CALLING(PEAK_CALLING.out.peaks)

    VARIANT_ANNOTATION(VARIANT_CALLING.out.vcf)

    VCF_POSTPROCESSING(VARIANT_ANNOTATION.out.vcf)

    emit:
    fastqc_html = QUALITY_CONTROL.out.fastqc_html
    maf_out     = VCF_POSTPROCESSING.out.maf
    vcf_out     = VCF_POSTPROCESSING.out.vcf
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

    VCF_POSTPROCESSING(VARIANT_ANNOTATION.out.vcf)

    emit:
    maf_out     = VCF_POSTPROCESSING.out.maf
    vcf_out     = VCF_POSTPROCESSING.out.vcf
}
