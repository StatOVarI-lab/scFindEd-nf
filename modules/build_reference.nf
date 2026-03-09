process BUILD_REFERENCE {

tag "${design_csv.baseName}"
  publishDir "${params.outdir}/reference_amplicon", mode: 'copy'

  conda "$projectDir/envs/reference.yml"

input:
    path design_csv

output:
  tuple path("ref_amplicons_sorted.fa"), path("ref_amplicons_sorted.fa.*")
  path "ref_amplicons_sorted.dict"



script:
"""
set -euo pipefail

Rscript $projectDir/bin/reference_from_designSummary.R ${design_csv}

# Sort
seqkit sort -N ref_amplicons.fa -o ref_amplicons_sorted.fa

# Index
bwa index -a bwtsw ref_amplicons_sorted.fa
samtools faidx ref_amplicons_sorted.fa
picard CreateSequenceDictionary \
  REFERENCE=ref_amplicons_sorted.fa \
  OUTPUT=ref_amplicons_sorted.dict
"""

}