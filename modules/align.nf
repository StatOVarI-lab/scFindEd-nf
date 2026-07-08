process ALIGN_BWA_MEM {

  tag "${sample_id}"
  publishDir "${params.outdir}/align", mode: 'copy'

  conda "$projectDir/envs/align.yml"

  cpus params.threads ?: 4

  input:
    tuple val(sample_id), path(r1), path(r2)
    tuple path(ref_fa), path(ref_idx)

  output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

  script:
"""
set -euo pipefail

bwa mem -t ${task.cpus} ${ref_fa} ${r1} ${r2} \\
  | samtools sort -@ ${task.cpus} -m 1G -o ${sample_id}.sorted.bam

samtools calmd -@ ${task.cpus} -b ${sample_id}.sorted.bam ${ref_fa} > ${sample_id}.bam

samtools index -@ ${task.cpus} ${sample_id}.bam
"""
}
