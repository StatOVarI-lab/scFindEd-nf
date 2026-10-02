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
    path("${sample_id}_samtools_stats.txt")
    path("${sample_id}_bamqc"), emit: plots

    script:
    """
    set -euo pipefail

    # 1) Basic alignment QC
    samtools flagstat \
        -@ ${task.cpus} \
        ${bam} \
        > ${sample_id}_flagstat.txt

    # 2) Reads mapped to each amplicon
    samtools idxstats \
        ${bam} \
        > ${sample_id}_idxstats.txt

    # 3) More detailed BAM statistics
    samtools stats \
        -@ ${task.cpus} \
        ${bam} \
        > ${sample_id}_samtools_stats.txt

    # 4) Generate BAM QC plots
    mkdir -p ${sample_id}_bamqc

    plot-bamstats \
        -p ${sample_id}_bamqc/ \
        ${sample_id}_samtools_stats.txt
    """
}