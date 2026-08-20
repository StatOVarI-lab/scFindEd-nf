#!/usr/bin/env Rscript

library(httr)
library(jsonlite)
library(stringr)
library(readr)

get_amplicon_sequences <- function(
  design_summary,
  identity_amplicons = character(0),
  verbose = FALSE
) {

  if (verbose) {
    message("Extracting amplicon sequences from Ensembl...")
  }

  # ---------------------------------------------------------
  # Column containing the amplicon ID
  # CHANGE THIS ONLY if your design CSV uses a different name
  # ---------------------------------------------------------
  amplicon_id_col <- "AmpID"

  if (!amplicon_id_col %in% names(design_summary)) {
    stop(
      "Amplicon ID column '",
      amplicon_id_col,
      "' was not found in design_summary.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------
  # Check that all requested identity amplicons actually exist
  # ---------------------------------------------------------
  if (length(identity_amplicons) > 0) {

    missing_identity_amplicons <- setdiff(
      identity_amplicons,
      design_summary[[amplicon_id_col]]
    )

    if (length(missing_identity_amplicons) > 0) {
      stop(
        "Identity amplicon(s) not found in design summary: ",
        paste(missing_identity_amplicons, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------
  # Select rows that should have sequences fetched from Ensembl.
  #
  # Identity amplicons are excluded BY ID.
  # Their position in the CSV does not matter.
  # ---------------------------------------------------------
  fetch_rows <- !design_summary[[amplicon_id_col]] %in% identity_amplicons

  temp_design_summary <- design_summary[fetch_rows, ]

  if (verbose) {
    if (length(identity_amplicons) > 0) {
      message(
        "Excluding identity amplicon(s): ",
        paste(identity_amplicons, collapse = ", ")
      )
    }

    message(
      "Fetching sequences for ",
      nrow(temp_design_summary),
      " amplicon(s)."
    )
  }

  # If everything was excluded, return the original design
  if (nrow(temp_design_summary) == 0) {
    return(design_summary)
  }

  # ---------------------------------------------------------
  # Check that only one reference genome is listed
  # among the amplicons being queried
  # ---------------------------------------------------------
  if (length(unique(temp_design_summary$Ref_Genome)) > 1) {

    stop(
      "Multiple reference genomes provided in design_summary.",
      call. = FALSE
    )

  } else {

    ref_genome <- temp_design_summary$Ref_Genome[1]

  }

  # ---------------------------------------------------------
  # Choose Ensembl server
  # ---------------------------------------------------------
  if (ref_genome %in% c("hg38", "GRCh38", "hg38_GR")) {

    server <- "https://rest.ensembl.org"

  } else if (ref_genome %in% c("hg19", "GRCh37", "hg19_GR")) {

    server <- "https://grch37.rest.ensembl.org"

  } else {

    stop(
      "Unsupported reference genome provided.",
      call. = FALSE
    )

  }

  ext <- "/sequence/region/human"

  # ---------------------------------------------------------
  # Build Ensembl regions
  # ---------------------------------------------------------
  regions_list <- paste0(
    stringr::str_sub(temp_design_summary$chr, 4, -1),
    ":",
    temp_design_summary$amplicon_start,
    "..",
    temp_design_summary$amplicon_end
  )

  body <- list(
    regions = regions_list
  )

  # ---------------------------------------------------------
  # Query Ensembl
  # ---------------------------------------------------------
  response <- httr::POST(
    url = paste0(server, ext),
    httr::content_type("application/json"),
    httr::accept("application/json"),
    body = jsonlite::toJSON(
      body,
      auto_unbox = TRUE
    )
  )

  httr::stop_for_status(response)

  results <- jsonlite::fromJSON(
    httr::content(
      response,
      as = "text",
      encoding = "UTF-8"
    )
  )

  output_sequences <- results$seq

  # Sanity check
  if (length(output_sequences) != nrow(temp_design_summary)) {
    stop(
      "Number of sequences returned by Ensembl does not match ",
      "the number of requested amplicons.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------
  # Preserve the ORIGINAL design row order.
  #
  # Only replace amplicon_seq for non-identity amplicons.
  # Identity rows remain exactly where they were.
  # ---------------------------------------------------------
  output_design <- design_summary

  if (!"amplicon_seq" %in% names(output_design)) {
    output_design$amplicon_seq <- NA_character_
  }

  output_design$amplicon_seq[fetch_rows] <- output_sequences

  if (verbose) {
    message(
      "Amplicon sequences extracted and appended to design_summary."
    )
  }

  return(output_design)
}


# ============================================================
# Command-line arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "get_amplicon_sequence.R",
      "<input_csv>",
      "<output_csv>",
      "[identity_amplicon1,identity_amplicon2,...]",
      "[--verbose]"
    ),
    call. = FALSE
  )
}

input_csv <- args[1]
output_csv <- args[2]

# Third argument is optional.
# Example:
# AMPL902654
#
# or:
# AMPL902654,AMPL902655
identity_amplicon_arg <- ""

if (length(args) >= 3 && args[3] != "--verbose") {
  identity_amplicon_arg <- args[3]
}

# Convert comma-separated string into an R vector
if (nzchar(identity_amplicon_arg)) {

  identity_amplicons <- strsplit(
    identity_amplicon_arg,
    ",",
    fixed = TRUE
  )[[1]]

  identity_amplicons <- trimws(identity_amplicons)

  identity_amplicons <- identity_amplicons[
    nzchar(identity_amplicons)
  ]

} else {

  identity_amplicons <- character(0)

}

verbose <- "--verbose" %in% args


# ============================================================
# Read design
# ============================================================

design_summary <- readr::read_csv(
  input_csv,
  show_col_types = FALSE
)


# ============================================================
# Get sequences
# ============================================================

design_summary_with_amplicon_seqs <- get_amplicon_sequences(
  design_summary = design_summary,
  identity_amplicons = identity_amplicons,
  verbose = verbose
)


# ============================================================
# Write output
# ============================================================

readr::write_csv(
  design_summary_with_amplicon_seqs,
  output_csv
)