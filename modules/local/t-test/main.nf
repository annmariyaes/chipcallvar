/*
=============================================================================
        MODULE: Perform t-test, permutation test on filtered variants 
                and find the significant genes/biomarkers
=============================================================================
*/


process TTEST {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/post_processing/dge_analysis", mode: 'copy'
    container "${params.R_CONTAINER}"
    // label 'process_RScript'
  
    input:
    tuple val(meta), path(csv)
    path(tpm)
    path(r_script)
    
    output:
    tuple val(meta), path("*.png"), emit: violin_plots, optional: true
    tuple val(meta), path("*.csv"), emit: csv
    
    script:
    """
    set -euo pipefail
    # Run the R script with proper arguments
    Rscript ${r_script} ${tpm} ${csv}
    """
}
