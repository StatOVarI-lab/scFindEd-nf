process GET_AMPLICON_SEQ {

    tag "${design_summary.simpleName}"

    conda "$projectDir/envs/get_amplicon_seq.yml"

    input:
    path design_summary
    val identity_amplicon

    output:
    path "${design_summary.simpleName}.with_amplicon_seq.csv"

    script:
    """
    Rscript ${projectDir}/bin/get_amplicon_seq.R \
        "${design_summary}" \
        "${design_summary.simpleName}.with_amplicon_seq.csv" \
        "${identity_amplicon}"
    """
}