process FASTQC {

    tag "${sample_id}"

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    conda "$projectDir/envs/fastqc.yml"

    cpus 2

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        ${r1} \
        ${r2}
    """
}