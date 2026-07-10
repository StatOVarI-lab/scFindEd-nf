process QC_PLOTS {

  tag "${sample_id}"
  publishDir "${params.outdir}/qc", mode: 'copy'

  conda "$projectDir/envs/qc.yml"

  cpus 2

  input:
    tuple val(sample_id), path(bam), path(bai)

  output:
    path("${sample_id}_flagstat.txt")
    path("${sample_id}_idxstats.txt")


  script:
  """
  set -euo pipefail

  # 1) Basic QC
  samtools flagstat -@ ${task.cpus} ${bam} > ${sample_id}_flagstat.txt
  samtools idxstats ${bam} > ${sample_id}_idxstats.txt
  """
}