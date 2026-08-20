process GET_AMPLICON_IDS {

    tag "${snp_info.baseName}"

    conda "${projectDir}/envs/get_amplicon_ids.yml"

    input:
    path snp_info
    path design_summary
    val identity_amplicon

    output:
    path "${snp_info.simpleName}.with_AmpIDs.csv"

    script:
    """
    Rscript ${projectDir}/bin/get_amplicon_ids.R \
        ${snp_info} \
        ${design_summary} \
        ${snp_info.simpleName}.with_AmpIDs.csv \
        "${identity_amplicon}" \
        --verbose
    """
}