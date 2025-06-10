nextflow.enable.dsl=2

include { CHIP_SEQ_FASTQ_VARIANT_CALLING } from './workflows/chipcallvar'
include { CHIP_SEQ_BAM_VARIANT_CALLING } from './workflows/chipcallvar'
include { CHIP_SEQ_VCF_VARIANT_ANNOTATION } from './workflows/chipcallvar'

/*
========================================================================================
   MAIN WORKFLOW - Variant Calling
========================================================================================
*/


workflow CHIPCALLVAR {

    // by default
    if (params.step == 'mapping') {
        ch_input = Channel
            .fromPath(params.SAMPLESHEET)
            .splitCsv(header: true)
            .map { row ->
                 def meta = [
                    id: "${row.sample}_${row.replicate}", 
                    patient: row.patient, 
                    sample: row.sample,
                    replicate: row.replicate.toInteger(),  
                    control: row.control
                ]
                // Handle both single-end and paired-end data
                def fastq_files = []
                
                // Check for different possible column names for FASTQ files
                if (row.containsKey('fastq_1') && row.fastq_1) {
                    fastq_files.add(file(row.fastq_1))

                     // Paired-end case with fastq_1 and fastq_2 columns
                    if (row.containsKey('fastq_2') && row.fastq_2) {
                        fastq_files.add(file(row.fastq_2))
                        meta.single_end = false
                    } else {
                        meta.single_end = true
                    }
                }                 
                return [ meta, fastq_files ]
            }
            // .view { it -> "$it" }
        CHIP_SEQ_FASTQ_VARIANT_CALLING(ch_input)
    }

    // 
    else if (params.step == "variant_calling") {
        ch_input = Channel
            .fromPath(params.SAMPLESHEET)
            .splitCsv(header: true)
            .map { row ->
                 def meta = [
                    id: "${row.sample}_${row.replicate}", 
                    patient: row.patient, 
                    sample: row.sample,
                    replicate: row.replicate.toInteger(),  
                    control: row.control,
                    single_end: row.single_end,
                ]         
                return [ [meta, file(row.bam), file(row.bai)] ]
            }
            // .view { it -> "$it" }
        CHIP_SEQ_BAM_VARIANT_CALLING(ch_input)
    }

    else if (params.step == "annotate") {
        ch_vcf = Channel
            .fromPath(params.SAMPLESHEET)
            .splitCsv(header: true)
            .filter { row -> row.vcf && file(row.vcf).exists() }
            .map { row ->
                def meta = [ id: row.sample ]
                [ meta, file(row.vcf) ]
            }
            // .view { it -> "$it" }
            
        CHIP_SEQ_VCF_VARIANT_ANNOTATION(ch_vcf)
    }
        
}


workflow {
    CHIPCALLVAR()
}


println """
\033[1;35m	ChIP-Seq Variant Calling Nextflow Pipeline\033[0m
\033[1;34m=============================================================================================================================\033[0m
\033[1;36m  Samplesheet       :\033[0m ${params.SAMPLESHEET}
\033[1;36m  Output Directory  :\033[0m ${params.OUTDIR}
\033[1;36m  Reference genome  :\033[0m ${params.REFERENCE_GENOME}
\033[1;36m  VEP Cache         :\033[0m ${params.VEP_CACHE}
\033[1;34m=============================================================================================================================\033[0m
"""



