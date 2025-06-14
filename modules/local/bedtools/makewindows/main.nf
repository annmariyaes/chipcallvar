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
    path "genome_chunks_*.bed", emit: intervals
    path "genome_chunks.bed", emit: bed_intervals

    script:
    """
    # bedtools makewindows -g ${fai} -w ${params.WINDOW_SIZE} > genome_chunks.bed
    awk -v FS="\t" -v OFS="\t" '{print \$1 FS "0" FS (\$2-1)}' ${fai} > genome_chunks.bed
    # Split using awk (more portable)
    awk '{print > "genome_chunks_" NR ".bed"}' genome_chunks.bed

    """
}
