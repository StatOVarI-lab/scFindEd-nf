#!/usr/bin/env Rscript

library(httr)
library(jsonlite)
library(stringr)
library(readr)

get_amplicon_sequences <- function(
  design_summary,
  verbose = FALSE
) {
  if (verbose) {
    message("Extracting amplicon sequences from Ensembl...")
  }

  temp_design_summary <- design_summary[-nrow(design_summary), ]

  # Check that only one reference genome is listed in design_summary
  if (length(unique(temp_design_summary$Ref_Genome)) > 1) {
    stop("Multiple reference genomes provided in design_summary.", call. = FALSE)
  } else {
    ref_genome <- temp_design_summary$Ref_Genome[1]
  }

  # Get correct server for the listed reference genome
  if (ref_genome %in% c("hg38", "GRCh38", "hg38_GR")) {
    server <- "https://rest.ensembl.org"
  } else if (ref_genome %in% c("hg19", "GRCh37", "hg19_GR")) {
    server <- "https://grch37.rest.ensembl.org"
  } else {
    stop("Unsupported reference genome provided.", call. = FALSE)
  }

  ext <- "/sequence/region/human"

  # Extract vector of regions from design_summary
  regions_list <- paste0(
    stringr::str_sub(temp_design_summary$chr, 4, -1), ":",
    temp_design_summary$amplicon_start, "..",
    temp_design_summary$amplicon_end
  )

  # Refactor as JSON for submission to Ensembl
  body <- list(regions = regions_list)

  # Send POST request
  response <- httr::POST(
    url = paste0(server, ext),
    httr::content_type("application/json"),
    httr::accept("application/json"),
    body = jsonlite::toJSON(body, auto_unbox = TRUE)
  )

  # Parse results
  httr::stop_for_status(response)

  results <- jsonlite::fromJSON(
    httr::content(response, as = "text", encoding = "UTF-8")
  )

  output_sequences <- results$seq

  # Overwrite design_summary$amplicon_seq with Ensembl-extracted sequences
  temp_design_summary$amplicon_seq <- output_sequences
  temp_design_summary[nrow(design_summary), ] <- design_summary[nrow(design_summary), ]

  if (verbose) {
    message("Amplicon sequences extracted and appended to design_summary.")
  }

  return(temp_design_summary)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: get_amplicon_sequence.R <input_csv> <output_csv> [--verbose]",
    call. = FALSE
  )
}

input_csv  <- args[1]
output_csv <- args[2]
verbose    <- "--verbose" %in% args

design_summary <- readr::read_csv(input_csv, show_col_types = FALSE)

design_summary_with_amplicon_seqs <- get_amplicon_sequences(
  design_summary = design_summary,
  verbose = verbose
)

readr::write_csv(design_summary_with_amplicon_seqs, output_csv)