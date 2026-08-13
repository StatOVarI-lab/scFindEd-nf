process REPORT_UTR_AMPLICON {

    tag "${sample_id}:${amp_id}"

    publishDir "${params.outdir}/utr_amplicon", mode: 'copy'

    cpus 1

    input:
    tuple val(sample_id), val(amp_id), path(amplicon_tsv), path(barcodes_readnames)

    output:
    tuple val(sample_id), val(amp_id),
        path("${sample_id}_${amp_id}_utr_reads.tsv.gz")

    script:
    """
    set -euo pipefail

    # Extract unique read names mapping to the UTR amplicon.
    # The amplicon TSV has a header and read_name is column 1.
    zcat ${amplicon_tsv} \
      | awk -F '\\t' 'NR > 1 { print \$1 }' \
      | LC_ALL=C sort -u \
      > utr_read_names.txt

    # Stream through the cell barcode/read-name file.
    # Only retain reads found in the UTR amplicon.
    {
      echo -e "sample_id\\tcell_barcode\\tamp_id\\tread_name"

      awk -F '\\t' \
        -v OFS='\\t' \
        -v sample="${sample_id}" \
        -v amp="${amp_id}" \
        '
        NR == FNR {
            wanted[\$1] = 1
            next
        }

        (\$1 in wanted) && \$2 != "" {
            print sample, \$2, amp, \$1
        }
        ' \
        utr_read_names.txt \
        <(zcat ${barcodes_readnames})

    } | gzip -c > ${sample_id}_${amp_id}_utr_reads.tsv.gz

    test -s ${sample_id}_${amp_id}_utr_reads.tsv.gz
    """
}