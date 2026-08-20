nextflow.enable.dsl = 2

include { GET_AMPLICON_SEQ }          from './modules/get_amplicon_seq'
include { GET_AMPLICON_IDS }          from './modules/get_amplicon_ids'
include { JOIN_DESIGN_AND_SNP }       from './modules/join_design_and_snp'
include { BUILD_REFERENCE }           from './modules/build_reference'
include { MAKE_AMPLICON_LIST }        from './modules/make_amplicon_list'
include { ALIGN_BWA_MEM }             from './modules/align'
include { QC_PLOTS }                  from './modules/qc_plots'
include { EXTRACT_AMPLICON_TSV }      from './modules/extract_amplicon_tsv'
include { EXTRACT_READNAMES }         from './modules/extract_readnames'
include { EXTRACT_READ_MUTATION }     from './modules/extract_read_mutation'
include { REPORT_identity_AMPLICON }  from './modules/report_identity_amplicon'
include { FASTQC }                    from './modules/fastqc'
include { MULTIQC }                   from './modules/multiqc'


workflow {

    /*
     * Required parameters
     */
    if (!params.design) {
        error "Please provide a design CSV file with --design"
    }

    if (!params.reads) {
        error "Please provide a glob pattern for reads with --reads"
    }

    if (!params.snp_info) {
        error "Please provide SNP info file with --snp_info"
    }


    /*
     * Identity amplicons
     *
     * Supports:
     *
     * --identity_amplicon AMPL902654
     *
     * or:
     *
     * --identity_amplicon AMPL902654,AMPL902655
     *
     * If no identity amplicon is supplied, this becomes an empty list.
     */
    identity_amplicons = params.identity_amplicon ?
        params.identity_amplicon
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it } :
        []


    /*
     * Keep the original comma-separated form for the R script.
     *
     * Example:
     *
     * "AMPL902654,AMPL902655"
     *
     * or "" if no identity amplicon was provided.
     */
    identity_amplicon_arg = params.identity_amplicon ?
        params.identity_amplicon.toString() :
        ''


    /*
     * Input channels
     */
    design_ch = Channel.fromPath(params.design)

    snp_info_ch = Channel.fromPath(params.snp_info)


    /*
     * Add amplicon sequences to the design CSV.
     *
     * Identity amplicons are passed to GET_AMPLICON_SEQ so that
     * their rows can be excluded by amplicon ID rather than by
     * assuming they are the last row(s) in the design file.
     */
    design_with_seq_ch = GET_AMPLICON_SEQ(
        design_ch,
        identity_amplicon_arg
    )


    /*
     * SNP branch
     */
    snp_with_ampid_ch = GET_AMPLICON_IDS(
        snp_info_ch,
        design_with_seq_ch,
        identity_amplicon_arg
    )

    joined_design_ch = JOIN_DESIGN_AND_SNP(
        design_with_seq_ch,
        snp_with_ampid_ch
    )


    /*
     * Build reference from the full design CSV
     */
    def (ref_bundle_ch, ref_dict_ch) = BUILD_REFERENCE(
        design_with_seq_ch
    )


    /*
     * Make list of SNP/HDR amplicons from the joined design
     */
    amplicon_list_file_ch = MAKE_AMPLICON_LIST(
        joined_design_ch
    )


    /*
     * Convert the amplicon list file into a channel of IDs
     */
    amplicon_ids_ch = amplicon_list_file_ch
        .splitText()
        .map { it.trim() }
        .filter { it }


    /*
     * Add optional identity amplicon(s) to the set of amplicons
     * used for EXTRACT_AMPLICON_TSV.
     *
     * Example:
     *
     * SNP amplicons:
     *   AMPL902647
     *   AMPL902649
     *
     * identity amplicons:
     *   AMPL902654
     *   AMPL902655
     *
     * resulting extraction set:
     *   AMPL902647
     *   AMPL902649
     *   AMPL902654
     *   AMPL902655
     */
    amplicon_ids_for_tsv_ch = identity_amplicons ?
        amplicon_ids_ch
            .mix(Channel.fromList(identity_amplicons))
            .unique() :
        amplicon_ids_ch


    /*
     * FASTQ input
     */
    reads_ch = Channel.fromFilePairs(
        params.reads,
        flat: true
    )


    FASTQC(reads_ch)

    /*
    * Collect FastQC results from all samples/read pairs
    * and generate one combined MultiQC report.
    */
    fastqc_files_ch = FASTQC.out.html
        .mix(FASTQC.out.zip)
        .collect()

    MULTIQC(fastqc_files_ch)


    /*
     * Extract R1 only for read-name / barcode extraction
     */
    r1_ch = reads_ch.map {
        sample_id,
        r1,
        r2 ->

        tuple(
            sample_id,
            r1
        )
    }


    /*
     * Extract read-name / cell-barcode mapping
     */
    readnames_ch = EXTRACT_READNAMES(
        r1_ch
    )


    /*
     * Align reads
     */
    bam_ch = ALIGN_BWA_MEM(
        reads_ch,
        ref_bundle_ch
    )


    /*
     * QC plots
     */
    QC_PLOTS(
        bam_ch
    )


    /*
     * Scatter over:
     *
     * SNP/HDR amplicons
     * +
     * optional identity amplicon(s)
     */
    scatter_ch = bam_ch
        .combine(amplicon_ids_for_tsv_ch)


    /*
     * Extract alignment TSV for every requested amplicon
     */
    amplicon_tsv_ch = EXTRACT_AMPLICON_TSV(
        scatter_ch
    )


    /*
     * SNP/HDR mutation branch
     *
     * Exclude ALL identity amplicons from mutation extraction.
     *
     * Their position in the design CSV is irrelevant.
     */
    mutation_tsv_ch = identity_amplicons ?
        amplicon_tsv_ch.filter { sample_id, amp_id, amplicon_tsv ->
            !identity_amplicons.contains(amp_id.toString())
        } :
        amplicon_tsv_ch


    /*
     * Combine each SNP/HDR amplicon TSV with the sample's
     * read-name / cell-barcode table
     */
    fragment_input_ch = mutation_tsv_ch
        .combine(
            readnames_ch,
            by: 0
        )
        .map {
            sample_id,
            amp_id,
            amplicon_tsv,
            barcodes_readnames ->

            tuple(
                sample_id,
                amp_id,
                amplicon_tsv,
                barcodes_readnames
            )
        }


    /*
     * Extract fragment mutation patterns
     */
    read_mutations_ch = EXTRACT_READ_MUTATION(
        fragment_input_ch
    )


    /*
     * Identity amplicon branch
     *
     * Select ALL requested identity amplicons.
     *
     * These are combined with the cell-aware read-name table,
     * but mutation extraction is not performed.
     */
    if (identity_amplicons) {

        identity_tsv_ch = amplicon_tsv_ch.filter {
            sample_id, amp_id, amplicon_tsv ->
                identity_amplicons.contains(amp_id.toString())
        }


        identity_input_ch = identity_tsv_ch
            .combine(
                readnames_ch,
                by: 0
            )
            .map {
                sample_id,
                amp_id,
                amplicon_tsv,
                barcodes_readnames ->

                tuple(
                    sample_id,
                    amp_id,
                    amplicon_tsv,
                    barcodes_readnames
                )
            }


        identity_report_ch = REPORT_identity_AMPLICON(
            identity_input_ch
        )
    }
}