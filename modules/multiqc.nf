process MULTIQC {

    tag "MultiQC"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    conda "$projectDir/envs/fastqc.yml"

    cpus 1

    input:
    path fastqc_files

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . \
        --filename multiqc_report.html
    """
}