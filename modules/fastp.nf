process FASTP {

    tag "${sample_id}"

    publishDir "${params.outdir}/qc/fastp", mode: 'copy'

    conda "$projectDir/envs/fastp.yml"

    cpus 4

    input:
    tuple val(sample_id),
          path(r1),
          path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fastq.gz"),
          path("${sample_id}_R2.trimmed.fastq.gz"),
          emit: reads

    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json

    script:
    """
    set -euo pipefail

    fastp \
        --in1 "${r1}" \
        --in2 "${r2}" \
        --out1 "${sample_id}_R1.trimmed.fastq.gz" \
        --out2 "${sample_id}_R2.trimmed.fastq.gz" \
        --thread ${task.cpus} \
        --cut_tail \
        --cut_tail_window_size 4 \
        --cut_tail_mean_quality 20 \
        --disable_adapter_trimming \
        --disable_quality_filtering \
        --disable_length_filtering \
        --disable_trim_poly_g \
        --html "${sample_id}.fastp.html" \
        --json "${sample_id}.fastp.json"
    """
}