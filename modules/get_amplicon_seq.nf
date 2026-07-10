process GET_AMPLICON_SEQ {
    tag "${design_summary.baseName}"
    conda "${projectDir}/envs/get_amplicon_seq.yml"
    // publishDir "${params.outdir}/design", mode: 'copy'

    input:
    path design_summary

    output:
    path "${design_summary.simpleName}.with_amplicon_seq.csv"

    script:
    """
    Rscript ${projectDir}/bin/get_amplicon_seq.R \
        ${design_summary} \
        ${design_summary.simpleName}.with_amplicon_seq.csv \
        --verbose
    """
}