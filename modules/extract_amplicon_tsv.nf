process EXTRACT_AMPLICON_TSV {

  tag "${sample_id}:${amp_id}"
  publishDir "${params.outdir}/amplicon_tsv/${sample_id}", mode: 'copy'

  conda "$projectDir/envs/qc.yml"   // needs samtools (and gzip). Reuse QC env or make a small one.
  cpus 1

  input:
    tuple val(sample_id), path(bam), path(bai), val(amp_id)

  output:
    tuple val(sample_id), val(amp_id), path("${sample_id}.${amp_id}.data.tsv.gz")

  script:
  """
  set -euo pipefail

  # Pull alignments mapping to this amplicon contig and extract columns:
  # QNAME, RNAME, POS, MAPQ, CIGAR, SEQ
  #
  # -F 4 excludes unmapped reads
  # If this amplicon has no reads, we still create an empty gz file.
  samtools view -F 4 ${bam} ${amp_id} \\
    | awk 'BEGIN {OFS="\\t"} {print \$1,\$3,\$4,\$5,\$6,\$10}' \\
    | gzip -c > ${sample_id}.${amp_id}.data.tsv.gz

  # Ensure file exists even if there were zero reads (samtools may exit 0 with no output anyway)
  test -f ${sample_id}.${amp_id}.data.tsv.gz
  """
}