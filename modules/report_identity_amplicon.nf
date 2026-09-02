process REPORT_identity_AMPLICON {

    tag "${sample_id}:${amp_id}"

    publishDir "${params.outdir}/identity_amplicon", mode: 'copy'

    conda "$projectDir/envs/report_identity_amplicon.yml"

    cpus 1

    input:
    tuple val(sample_id),
          val(amp_id),
          path(amplicon_tsv),
          path(barcodes_readnames)

    output:
    tuple val(sample_id),
          val(amp_id),
          path("${sample_id}_${amp_id}_identity_reads.tsv.gz")

    script:
    """
    set -euo pipefail

    Rscript ${projectDir}/bin/report_identity_amplicon.R \
        "${sample_id}" \
        "${amp_id}" \
        "${amplicon_tsv}" \
        "${barcodes_readnames}" \
        "${sample_id}_${amp_id}_identity_reads.tsv"

    gzip -f "${sample_id}_${amp_id}_identity_reads.tsv"

    test -s "${sample_id}_${amp_id}_identity_reads.tsv.gz"
    """
}