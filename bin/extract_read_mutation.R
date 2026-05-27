#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    "Usage: group_fragment_patterns.R <sample_id> <amp_id> <amplicon_tsv.gz> <barcodes_readnames.tsv.gz> <output_tsv>",
    call. = FALSE
  )
}

sample_id <- args[1]
amp_id_arg <- args[2]
input_file <- args[3]
barcode_file <- args[4]
output_file <- args[5]

amp <- readr::read_tsv(input_file, show_col_types = FALSE)

barcodes <- readr::read_tsv(
  barcode_file,
  col_names = c("read_name", "cell_barcode"),
  show_col_types = FALSE
)

amp <- amp %>%
  left_join(barcodes, by = "read_name") %>%
  filter(!is.na(cell_barcode), cell_barcode != "")

empty_out <- tibble(
  sample_id = character(),
  cell_barcode = character(),
  amp_id = character(),
  r1_pos = integer(),
  r1_mapq = integer(),
  r1_cigar = character(),
  r1_flag = integer(),
  r1_nm = integer(),
  r1_md = character(),
  r1_seq = character(),
  r2_pos = integer(),
  r2_mapq = integer(),
  r2_cigar = character(),
  r2_flag = integer(),
  r2_nm = integer(),
  r2_md = character(),
  r2_seq = character(),
  n_fragments = integer()
)

if (nrow(amp) == 0) {
  readr::write_tsv(empty_out, output_file)
  quit(save = "no")
}

amp_grouped_input <- amp %>%
  mutate(
    flag = as.integer(flag),
    nm = as.integer(nm),

    is_paired = bitwAnd(flag, 1) != 0,
    is_proper_pair = bitwAnd(flag, 2) != 0,
    is_read1 = bitwAnd(flag, 64) != 0,
    is_read2 = bitwAnd(flag, 128) != 0,

    is_secondary = bitwAnd(flag, 256) != 0,
    is_supplementary = bitwAnd(flag, 2048) != 0
  ) %>%
  filter(
    is_paired,
    is_proper_pair,
    !is_secondary,
    !is_supplementary
  )

if (nrow(amp_grouped_input) == 0) {
  readr::write_tsv(empty_out, output_file)
  quit(save = "no")
}

r1 <- amp_grouped_input %>%
  filter(is_read1) %>%
  select(
    read_name,
    cell_barcode,
    amp_id,
    r1_pos = pos,
    r1_mapq = mapq,
    r1_cigar = cigar,
    r1_flag = flag,
    r1_nm = nm,
    r1_md = md,
    r1_seq = seq
  )

r2 <- amp_grouped_input %>%
  filter(is_read2) %>%
  select(
    read_name,
    cell_barcode,
    amp_id,
    r2_pos = pos,
    r2_mapq = mapq,
    r2_cigar = cigar,
    r2_flag = flag,
    r2_nm = nm,
    r2_md = md,
    r2_seq = seq
  )

paired <- inner_join(
  r1,
  r2,
  by = c("read_name", "cell_barcode", "amp_id")
) %>%
  mutate(
    sample_id = sample_id
  )

if (nrow(paired) == 0) {
  readr::write_tsv(empty_out, output_file)
  quit(save = "no")
}

out <- paired %>%
  group_by(
    sample_id,
    cell_barcode,
    amp_id,
    r1_pos,
    r1_mapq,
    r1_cigar,
    r1_flag,
    r1_nm,
    r1_md,
    r1_seq,
    r2_pos,
    r2_mapq,
    r2_cigar,
    r2_flag,
    r2_nm,
    r2_md,
    r2_seq
  ) %>%
  summarise(
    n_fragments = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_fragments))

readr::write_tsv(out, output_file)