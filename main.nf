nextflow.enable.dsl=2

include { BUILD_REFERENCE } from './modules/build_reference'
include { MAKE_AMPLICON_LIST } from './modules/make_amplicon_list'
include { ALIGN_BWA_MEM } from './modules/align'
include { QC_PLOTS } from './modules/qc_plots'
include { EXTRACT_AMPLICON_TSV } from './modules/extract_amplicon_tsv'
include { EXTRACT_READNAMES } from './modules/extract_readnames'


workflow {
    if (! params.design) error "Please provide a design CSV file with --design"
    if (! params.reads) error "Please provide a glob pattern for reads with --reads"

 //Creates a channel that emits the design CSV file path
 design_ch = Channel.fromPath(params.design)
 //Runs BUILD_REFERENCE on the design CSV and captures its two output channels:
 def (ref_bundle_ch, ref_dict_ch) = BUILD_REFERENCE(design_ch)
 // MAKE_AMPLICON_LIST -> channel of amplicon IDs
 amplicon_list_file_ch = MAKE_AMPLICON_LIST(design_ch)
 //amplicon_list_file_ch.view { "LIST_FILE: $it" }

 amplicon_ids_ch = amplicon_list_file_ch
  .splitText()
  .map { it.trim() }
  .filter { it }

 //amplicon_ids_ch.view { "AMP_ID: '${it}'" }
 //Scans for read pairs matching your glob and emits tuples like:
 reads_ch = Channel.fromFilePairs(params.reads, flat: true)
 //readnames from R1
 r1_ch = reads_ch.map { sample_id, r1, r2 ->
        tuple(sample_id, r1)
    }

  readnames_ch = EXTRACT_READNAMES(r1_ch)
 //Aligns each read pair to the reference and emits
 bam_ch = ALIGN_BWA_MEM(reads_ch, ref_bundle_ch)
 //bam_ch.view { "BAM_TUPLE: $it" }
 QC_PLOTS(bam_ch)
 // Scatter: create one item per (sample, amplicon)
 scatter_ch = bam_ch
  .combine(amplicon_ids_ch)
  
 //scatter_ch.view { "SCATTER: $it" }

 amplicon_tsv_ch = EXTRACT_AMPLICON_TSV(scatter_ch)
     
 }
