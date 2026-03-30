#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(cli)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: join_design_and_snp_info.R <design_summary_csv> <snp_info_csv> <output_csv>",
    call. = FALSE
  )
}

design_file <- args[1]
snp_file <- args[2]
output_file <- args[3]

design_summary <- readr::read_csv(design_file, show_col_types = FALSE)
snp_info <- readr::read_csv(snp_file, show_col_types = FALSE)

# Check required columns
if (!"AmpID" %in% colnames(design_summary)) {
  cli::cli_abort("Column 'AmpID' not found in design_summary.")
}

if (!"AmpID" %in% colnames(snp_info)) {
  cli::cli_abort("Column 'AmpID' not found in SNP_info.")
}

# Clean AmpID
design_summary <- design_summary %>%
  mutate(AmpID = trimws(as.character(AmpID)))

snp_info <- snp_info %>%
  mutate(AmpID = trimws(as.character(AmpID)))

# Identify missing AmpIDs
missing_rows <- snp_info %>%
  filter(is.na(AmpID) | AmpID == "")

if (nrow(missing_rows) > 0) {
  warning_msg <- paste(
    apply(missing_rows, 1, function(row) {
      paste0(row["Gene"], " - ", row["SNP_number"])
    }),
    collapse = "\n"
  )

  cli::cli_warn(c(
    "Some SNPs could not be matched to an AmpID and will be excluded:",
    "!" = warning_msg
  ))
}

# Keep only valid SNPs
snp_info_filtered <- snp_info %>%
  filter(!is.na(AmpID) & AmpID != "")

# Check that remaining AmpIDs exist in design_summary
missing_in_design <- setdiff(unique(snp_info_filtered$AmpID), unique(design_summary$AmpID))
if (length(missing_in_design) > 0) {
  cli::cli_abort(c(
    "Some SNP AmpIDs were not found in design_summary:",
    "x" = paste(missing_in_design, collapse = ", ")
  ))
}

# Join
joined <- dplyr::inner_join(design_summary, snp_info_filtered, by = "AmpID")

readr::write_csv(joined, output_file)