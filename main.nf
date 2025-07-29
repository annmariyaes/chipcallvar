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
                    unique_id: "${row.sample}_${row.replicate}", 
                    id: row.id, 
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
                    unique_id: "${row.sample}_${row.replicate}", 
                    id: row.id, 
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
                def meta = [ 
                             sample: row.sample, 
                             id: row.id,
                             caller: row.caller ]
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

\033[1;36m  Starting step     :\033[0m ${params.STEP}
\033[1;36m  Variant callers   :\033[0m ${params.tools}
\033[1;36m  Skipped tools     :\033[0m ${params.skip_tools}

\033[1;34m  Reference genome:
\033[1;36m  Genome Assembly   :\033[0m ${params.ASSEMBLY}
\033[1;36m  Fasta             :\033[0m ${params.REFERENCE_GENOME}
\033[1;36m  Fasta index       :\033[0m ${params.GENOME_FAI}

\033[1;34m  Variant annotation database:
\033[1;36m  VEP Cache         :\033[0m ${params.VEP_CACHE}

\033[1;34m  Population databases:
\033[1;36m  Panel of Normals  :\033[0m ${params.PON}
\033[1;36m  Germline database :\033[0m ${params.GNOMAD}
\033[1;36m  SNVs database   :\033[0m ${params.DBSNP}
\033[1;36m  INDELs database   :\033[0m ${params.MILLS_1000G}    

\033[1;34m  Variant Filtering params:
\033[1;36m  Depth            :\033[0m ${params.DP}
\033[1;36m  AF (1000 genome) :\033[0m ${params.AF}
\033[1;36m  AF (gnomAD)      :\033[0m ${params.gnomADe_AF}
\033[1;36m  VAF              :\033[0m ${params.VAF}
\033[1;34m=============================================================================================================================\033[0m
"""



