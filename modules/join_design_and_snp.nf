process JOIN_DESIGN_AND_SNP {
    tag "${design_summary.baseName}"
    conda "$projectDir/envs/reference.yml"
    publishDir "${params.outdir}/design", mode: 'copy'

    input:
    path design_summary
    path snp_info

    output:
    path "${design_summary.simpleName}.joined.csv"

    script:
    """
    join_design_and_snp.R \
        ${design_summary} \
        ${snp_info} \
        ${design_summary.simpleName}.joined.csv
    """
}