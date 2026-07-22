#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]
compress <- args[3] == "TRUE"


extract_readnames <- function(input, output, compress) {
  con <- if (grepl("\\.gz$", input)) gzfile(input, "rt") else file(input, "r")

  out <- if (compress) {
    gzfile(output, "wt")
  } else {
    file(output, "wt")
  }

  on.exit({
    close(con)
    close(out)
  }, add = TRUE)

  i <- 0L
  repeat {
    lines <- readLines(con, n = 400000, warn = FALSE)
    if (length(lines) == 0) break

    idx <- seq(1, length(lines), by = 4)
    headers <- lines[idx]

    parsed <- sub("^@", "", headers)
    readname <- sub(" RG:Z:.*$", "", parsed)
    barcode <- ifelse(grepl(" RG:Z:", headers),
                      sub("^.* RG:Z:", "", headers),
                      "")

    writeLines(paste(readname, barcode, sep = "\t"), out)

    i <- i + length(idx)
    message("Processed reads: ", i)
  }
}

extract_readnames(input, output, compress)