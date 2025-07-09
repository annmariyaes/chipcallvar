/*
=============================================================================
        SUBWORKFLOW: 
   - Convert VCF files to MAF files
   - Plot MAF files for cohort analysis
   - Filter and perform t-test and permutation to find significant genes
     (RNA-seq integration)   
=============================================================================
*/

include { VCF2MAF } from '../../../../modules/local/vcf2maf'
include { MAFTOOLS } from '../../../../modules/local/maftools'
include { TTEST } from '../../../../modules/local/t-test'

workflow DOWNSTREAM_ANALYSIS {
    take:
    ch_vcf // channel: [ meta, vcf ]
    ch_tpm // From RNA-seq
    
    main:
    vcf2maf = VCF2MAF(ch_vcf)
    ch_maf = vcf2maf.maf.map { meta, maf ->
        [meta, maf.parent]
    }
    // .view { it -> "$it" }
    r_script_maftools = file("${projectDir}/bin/maftools_plotting.R")
    MAFTOOLS(ch_maf, r_script_maftools)
    
    // give option to perform differential expression analysis (t-tests) or not
    if (!(params.skip_tools && params.skip_tools.split(',').contains('dge'))) {
        r_script_dge = file("${projectDir}/bin/differential_expression_analysis.R")
        TTEST(MAFTOOLS.out.csv, ch_tpm, r_script_dge)
        
        ch_violin_plots = TTEST.out.violin_plots
        ch_ttest_csv = TTEST.out.csv
    } 
    else {
        ch_violin_plots = Channel.empty()
        ch_ttest_csv = Channel.empty()
    }
    
    emit:
    maf   = vcf2maf.maf  // [ meta, maf ]
    plots = MAFTOOLS.out.plots // waterfall and summary plots
    variants_csv = MAFTOOLS.out.csv // filtered, recurrent variants
    violin_plots = ch_violin_plots // violin plot of significant genes
    csv = ch_ttest_csv // filtered significant variants info
}
