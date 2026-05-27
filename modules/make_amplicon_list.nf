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
  if (!("amplicon_seq" %in% names(df))) stop("Column amplicon_seq not found")

  df\$AmpID <- trimws(df\$AmpID)
  df[["amplicon_seq"]] <- trimws(df[["amplicon_seq"]])

  df <- df[
    !is.na(df\$AmpID) &
    df\$AmpID != "" &
    !is.na(df[["amplicon_seq"]]) &
    df[["amplicon_seq"]] != "" &
    df[["amplicon_seq"]] != "NA",
  ]

  amps <- unique(df\$AmpID)
  writeLines(amps, "amplicon_list.txt")
'
  """
}