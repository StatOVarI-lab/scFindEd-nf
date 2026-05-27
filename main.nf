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
    if (! params.design) error "Please provide a design CSV file with --design"
    if (! params.reads) error "Please provide a glob pattern for reads with --reads"
    if ( !params.snp_info ) error "Please provide SNP info file with --snp_info"

    design_ch = Channel.fromPath(params.design)
    snp_info_ch       = Channel.fromPath(params.snp_info)

    // Add amplicon sequences to the design CSV, producing a new CSV with an added column.
    design_with_seq_ch = GET_AMPLICON_SEQ(design_ch)

    // SNP branch
    snp_with_ampid_ch = GET_AMPLICON_IDS(snp_info_ch, design_with_seq_ch)
    JOIN_DESIGN_AND_SNP(design_with_seq_ch, snp_with_ampid_ch)

    //Runs BUILD_REFERENCE on the design CSV and captures its two output channels:
    def (ref_bundle_ch, ref_dict_ch) = BUILD_REFERENCE(design_with_seq_ch)
    // MAKE_AMPLICON_LIST -> channel of amplicon IDs
    amplicon_list_file_ch = MAKE_AMPLICON_LIST(design_with_seq_ch)
    //amplicon_list_file_ch.view { "LIST_FILE: $it" }

    amplicon_ids_ch = amplicon_list_file_ch
     .splitText()
     .map { it.trim() }
     .filter { it }

    //Scans for read pairs matching your glob and emits tuples like:
    reads_ch = Channel.fromFilePairs(params.reads, flat: true)
    //readnames from R1
    r1_ch = reads_ch.map { sample_id, r1, r2 ->
        tuple(sample_id, r1)
    }

    readnames_ch = EXTRACT_READNAMES(r1_ch)
    //Aligns each read pair to the reference and emits
    bam_ch = ALIGN_BWA_MEM(reads_ch, ref_bundle_ch)

    QC_PLOTS(bam_ch)
    // Scatter: create one item per (sample, amplicon)
    scatter_ch = bam_ch
     .combine(amplicon_ids_ch)

    amplicon_tsv_ch = EXTRACT_AMPLICON_TSV(scatter_ch)

    fragment_input_ch = amplicon_tsv_ch
     .combine(readnames_ch, by: 0)
     .map { sample_id, amp_id, amplicon_tsv, barcodes_readnames ->
     tuple(sample_id, amp_id, amplicon_tsv, barcodes_readnames)
  }

    read_mutations_ch = EXTRACT_READ_MUTATION(fragment_input_ch)

     
 }