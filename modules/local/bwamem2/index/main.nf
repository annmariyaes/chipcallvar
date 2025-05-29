//chipcallvar/modules/nf-core/bwa/index/main.nf

nextflow.enable.dsl=2

process BWAMEM2_INDEX {
    tag { new File(params.REFERENCE_GENOME).getName() }
    publishDir "${params.OUTDIR}/reference", mode: 'copy'    
    container "${params.BWAMEM2_CONTAINER}"    
    label 'process_high'
    
    input:
    path reference

    output:
    tuple path(reference), path("${reference}.amb"), path("${reference}.ann"), path("${reference}.bwt"), path("${reference}.pac"), path("${reference}.sa"), emit: index

    script:
    """
    bwa-mem2 index ${reference}
    """
}
