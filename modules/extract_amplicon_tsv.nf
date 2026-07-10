process EXTRACT_AMPLICON_TSV {

  tag "${sample_id}:${amp_id}"
  // publishDir "${params.outdir}/amplicon_tsv/${sample_id}", mode: 'copy'

  conda "$projectDir/envs/qc.yml"
  cpus 1

  input:
    tuple val(sample_id), path(bam), path(bai), val(amp_id)

  output:
    tuple val(sample_id), val(amp_id), path("${sample_id}.${amp_id}.data.tsv.gz")

  script:
"""
set -euo pipefail

# Pull alignments mapping to this amplicon contig and extract:
# QNAME, RNAME, POS, MAPQ, CIGAR, FLAG, NM, MD, SEQ, QUAL
#
# -F 4 excludes unmapped reads.
# NM and MD are expected to be present because ALIGN_BWA_MEM now runs samtools calmd.

{
  echo -e "read_name\\tamp_id\\tpos\\tmapq\\tcigar\\tflag\\tnm\\tmd\\tseq\\tqual"

  samtools view -F 4 ${bam} ${amp_id} \\
    | awk 'BEGIN {OFS="\\t"}
      {
        nm="NA";
        md="NA";

        for (i=12; i<=NF; i++) {
          if (\$i ~ /^NM:i:/) {
            nm=\$i;
            sub(/^NM:i:/, "", nm);
          }

          if (\$i ~ /^MD:Z:/) {
            md=\$i;
            sub(/^MD:Z:/, "", md);
          }
        }

        print \$1, \$3, \$4, \$5, \$6, \$2, nm, md, \$10, \$11
      }'
} | gzip -c > ${sample_id}.${amp_id}.data.tsv.gz

test -f ${sample_id}.${amp_id}.data.tsv.gz
"""
}