#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(Rcpp)
})

# ---- Command-line options ----
option_list <- list(
  make_option("--input", type="character", help="Input AIRR TSV file"),
  make_option("--output", type="character", help="Output representative clone file"),
  make_option("--log", type="character", help="Text log file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---- Load Rcpp function to calculate mutations ----
Rcpp::sourceCpp(file = "/mnt/bin/rscripts/clones/allele_diff.cpp")

# ---- Read input AIRR TSV file ----
data <- read_tsv(opt$input, show_col_types = FALSE)

pre_selection <- nrow(data)

# ---- Calculate mutations ----
data$mut <- allele_diff(data$sequence_alignment, data$germline_alignment_d_mask, parallel = TRUE, return_count = TRUE)

# data$mut <- sapply(seq_len(nrow(data)), function(j) {
#   allele_diff(c(data$sequence_alignment[j], data$germline_alignment_d_mask[j]))
# })

# ---- Select representative per clone (fewest mutations) ----
data <- data %>%
  group_by(clone_id) %>%
  mutate(clone_size = n()) %>%
  slice_min(order_by = mut, with_ties = FALSE)

# ---- Write output files ----
write_tsv(data, file = opt$output)

if(!is.null(opt$log)){
    log_lines <- c(
    "START> Single Clone representative selection",
    sprintf("  INPUT> %s", opt$input),
    sprintf("  PRE-SELECTION> %s", pre_selection),
    sprintf("  PASS> %s", nrow(data)),
    "  END> Single Clone representative selection",
    ""
    )

    writeLines(log_lines, con = opt$log)

}
