
nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp.nf'
include { CHECK_STRANDNESS } from './modules/check_strandness.nf'
include { HISAT2_INDEX_REFERENCE ; HISAT2_ALIGN ; EXTRACT_SPLICE_SITES ; EXTRACT_EXONS } from './modules/hisat2.nf'
include { SAMTOOLS ; SAMTOOLS_MERGE } from './modules/samtools.nf'
include { CUFFLINKS } from './modules/cufflinks.nf'

log.info """\
         RNAseq analysis using NextFlow 
         =============================
         genome: ${params.reference_genome}
         annot : ${params.reference_annotation}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()
 
params.outdir = 'results'

workflow {

    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
    CHECK_STRANDNESS( read_pairs_ch, params.reference_cdna, params.reference_annotation_ensembl )
    FASTP( read_pairs_ch )
    EXTRACT_EXONS( params.reference_annotation )
    EXTRACT_SPLICE_SITES( params.reference_annotation )
    HISAT2_INDEX_REFERENCE( params.reference_genome, EXTRACT_EXONS.out, EXTRACT_SPLICE_SITES.out )
    HISAT2_ALIGN( read_pairs_ch, HISAT2_INDEX_REFERENCE.out, CHECK_STRANDNESS.out.first() )
    SAMTOOLS( HISAT2_ALIGN.out.sample_sam )
    CUFFLINKS( CHECK_STRANDNESS.out, SAMTOOLS.out.sample_bam, params.reference_annotation )
}
