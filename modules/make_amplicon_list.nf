process MAKE_AMPLICON_LIST {

  tag "${design_csv.baseName}"
  publishDir "${params.outdir}/amplicons", mode: 'copy'
  conda "$projectDir/envs/reference.yml"

  input:
    path design_csv

  output:
    path "amplicon_list.txt"

  script:
  """
  set -euo pipefail

  Rscript -e '
    df <- read.csv("${design_csv}", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

    if (!("AmpID" %in% names(df))) stop("Column AmpID not found")
    if (!("Amplicon sequence" %in% names(df))) stop("Column Amplicon sequence not found")

    df\$AmpID <- trimws(df\$AmpID)
    df[["Amplicon sequence"]] <- trimws(df[["Amplicon sequence"]])

    df <- df[
      !is.na(df\$AmpID) &
      df\$AmpID != "" &
      !is.na(df[["Amplicon sequence"]]) &
      df[["Amplicon sequence"]] != "" &
      df[["Amplicon sequence"]] != "NA",
    ]

    amps <- unique(df\$AmpID)
    writeLines(amps, "amplicon_list.txt")
  '
  """
}