nextflow.enable.dsl=2

include { GET_AMPLICON_SEQ }    from './modules/get_amplicon_seq'
include { GET_AMPLICON_IDS }    from './modules/get_amplicon_ids'
include { JOIN_DESIGN_AND_SNP } from './modules/join_design_and_snp'
include { BUILD_REFERENCE } from './modules/build_reference'
include { MAKE_AMPLICON_LIST } from './modules/make_amplicon_list'
include { ALIGN_BWA_MEM } from './modules/align'
include { QC_PLOTS } from './modules/qc_plots'
include { EXTRACT_AMPLICON_TSV } from './modules/extract_amplicon_tsv'
include { EXTRACT_READNAMES } from './modules/extract_readnames'
include { EXTRACT_READ_MUTATION } from './modules/extract_read_mutation'


workflow {

    design_ch = Channel.fromPath(params.design, checkIfExists: true)
    snp_info_ch = Channel.fromPath(params.snp_info, checkIfExists: true)

    // Add amplicon sequences to the design CSV, producing a new CSV with an added column.
    design_with_seq_ch = GET_AMPLICON_SEQ(design_ch)

  // SNP branch
    snp_with_ampid_ch = GET_AMPLICON_IDS(snp_info_ch, design_with_seq_ch)
    joined_design_ch = JOIN_DESIGN_AND_SNP(design_with_seq_ch, snp_with_ampid_ch)
  // Runs BUILD_REFERENCE on the full design CSV
    def (ref_bundle_ch, ref_dict_ch) = BUILD_REFERENCE(design_with_seq_ch)

  // MAKE_AMPLICON_LIST now uses only joined SNP/design amplicons
   amplicon_list_file_ch = MAKE_AMPLICON_LIST(joined_design_ch)

   amplicon_ids_ch = amplicon_list_file_ch
     .splitText()
     .map { it.trim() }
     .filter { it }

  // Add optional UTR barcode amplicon to amplicon TSV extraction only
  amplicon_ids_for_tsv_ch = params.utr_amplicon ?
  amplicon_ids_ch.mix(Channel.of(params.utr_amplicon)).unique() :
  amplicon_ids_ch

  reads_ch = Channel.fromFilePairs(params.reads, flat: true, checkIfExists: true)
  r1_ch = reads_ch.map { sample_id, r1, r2 ->
  tuple(sample_id, r1)
 }

  readnames_ch = EXTRACT_READNAMES(r1_ch)

  bam_ch = ALIGN_BWA_MEM(reads_ch, ref_bundle_ch.first())
  QC_PLOTS(bam_ch)

  // Scatter over SNP/HDR amplicons + optional UTR amplicon
  scatter_ch = bam_ch
     .combine(amplicon_ids_for_tsv_ch)

  amplicon_tsv_ch = EXTRACT_AMPLICON_TSV(scatter_ch)

  // Exclude UTR amplicon from fragment mutation extraction
  mutation_tsv_ch = params.utr_amplicon ?
  amplicon_tsv_ch.filter { sample_id, amp_id, amplicon_tsv ->
    amp_id != params.utr_amplicon
  } :
  amplicon_tsv_ch

  fragment_input_ch = mutation_tsv_ch
     .combine(readnames_ch, by: 0)
     .map { sample_id, amp_id, amplicon_tsv, barcodes_readnames ->
    tuple(sample_id, amp_id, amplicon_tsv, barcodes_readnames)
  }

read_mutations_ch = EXTRACT_READ_MUTATION(fragment_input_ch)
 }