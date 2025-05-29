include { FASTQC } from '../../../../modules/local/fastqc'


/*
=============================================================================
	SUBWORKFLOW: Quality Control
=============================================================================
*/

workflow QUALITY_CONTROL {
    take:
    ch_input // channel: [ meta, [ reads ] ]
    
    main:
    FASTQC(ch_input)
    
    ch_multiqc_files = FASTQC.out.zip.map{ meta, zip -> zip }.collect()
    
    emit:
    fastqc_html     = FASTQC.out.html
    fastqc_zip      = FASTQC.out.zip
}
