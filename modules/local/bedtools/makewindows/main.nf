/*
=============================================================================
        MODULE: BEDTOOLS makewindows
To speed up the variant calling processes, the reference is chopped into smaller pieces. 
The variant calling is done by this intervals, and the different resulting VCFs are then merged. 
This can parallelize the variant calling processes, and push down the variant calling wall clock time significantly.
=============================================================================
*/


process CREATE_INTERVALS_BED {
    publishDir "${params.OUTDIR}/preprocessing/intervals", mode: 'copy'
    container "${params.BEDTOOLS_CONTAINER}"

    input:
    path fai

    output:
    path "genome_${params.WINDOW_SIZE}bp_intervals.bed", emit: intervals

    script:
    """
    bedtools makewindows -g ${fai} -w ${params.WINDOW_SIZE} > genome_${params.WINDOW_SIZE}bp_intervals.bed
    """
}
