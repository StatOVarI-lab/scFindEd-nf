process REPORT_UTR_AMPLICON {

  tag "${sample_id}:${amp_id}"

  publishDir "${params.outdir}/utr_amplicon",
    mode: 'copy',
    pattern: "*.tsv.gz"

  conda "$projectDir/envs/qc.yml"

  input:
  tuple val(sample_id),
        val(amp_id),
        path(amplicon_tsv),
        path(barcodes_readnames)

  output:
  tuple val(sample_id),
        val(amp_id),
        path("${sample_id}_${amp_id}_utr_reads.tsv.gz")

  script:
  """
  set -euo pipefail

  Rscript ${projectDir}/bin/report_utr_amplicon.R \
    "${sample_id}" \
    "${amp_id}" \
    "${amplicon_tsv}" \
    "${barcodes_readnames}" \
    "${sample_id}_${amp_id}_utr_reads.tsv"

  gzip -f "${sample_id}_${amp_id}_utr_reads.tsv"
  """
}