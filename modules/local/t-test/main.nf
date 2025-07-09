/*
=============================================================================
        MODULE: Perform t-test, permutation test on filtered variants 
                and find the significant genes/biomarkers
=============================================================================
*/


process TTEST {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/analysis", mode: 'copy'
    container "${params.R_CONTAINER}"
    
    input:
    tuple val(meta), path(csv)
    tuple val(meta), path(tpm)
    path(r_script)
    
    output:
    tuple val(meta), path("*.png"), emit: violin_plots
    tuple val(meta), path("*.csv"), emit: csv
    
    script:
    """
    set -euo pipefail
    # Run the R script with proper arguments
    Rscript ${r_script}
    """
}
