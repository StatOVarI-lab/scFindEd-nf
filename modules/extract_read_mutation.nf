process EXTRACT_READ_MUTATION {

  tag "${sample_id}:${amp_id}"
  publishDir "${params.outdir}/fragment_mutations", mode: 'copy'

  conda "$projectDir/envs/read_mutations.yml"

  cpus 1

  input:
  tuple val(sample_id), val(amp_id), path(amplicon_tsv), path(barcodes_readnames)

  output:
    tuple val(sample_id), val(amp_id), path("${sample_id}_${amp_id}_fragment_mutations.tsv.gz")

  script:
  """
  set -euo pipefail

Rscript ${projectDir}/bin/extract_read_mutation.R \
  ${sample_id} \
  ${amp_id} \
  ${amplicon_tsv} \
  ${barcodes_readnames} \
  ${sample_id}_${amp_id}_fragment_mutations.tsv

  gzip -f ${sample_id}_${amp_id}_fragment_mutations.tsv
  """
}