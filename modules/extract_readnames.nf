process EXTRACT_READNAMES {

    // publishDir "${params.outdir}/readnames", mode: 'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(input_fastq)

    output:
    tuple val(sample_id), path("${sample_id}.readnames.tsv${params.compress ? '.gz' : ''}")

    script:
    def outfile = "${sample_id}.readnames.tsv" + (params.compress ? ".gz" : "")
    def compress_flag = params.compress ? "TRUE" : "FALSE"

    """
    Rscript ${projectDir}/bin/extract_readnames.R \
        "${input_fastq}" \
        "${outfile}" \
        "${compress_flag}"

    # Fail loudly on a truncated stream rather than publishing a short file that
    # decompresses partway and silently drops reads downstream.
    if [ "${compress_flag}" = "TRUE" ]; then
        gzip -t "${outfile}"
    fi
    """
}