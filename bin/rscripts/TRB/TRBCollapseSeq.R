#!/usr/bin/env Rscript

# Load required libraries
suppressMessages(library(piglet))
suppressMessages(library(tigger))
suppressMessages(library(dplyr))
library(stringr)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("-m", "--makedb_file"), type="character", help="Path to the makedb file (input TSV)", metavar="FILE"),
  make_option(c("-v", "--v_germline"), type="character", help="Path to the V germline FASTA file", metavar="FILE"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$makedb_file) || is.null(opt$v_germline)) {
  print_help(opt_parser)
  quit(status=1)
}

# Assign variables
makedb_file <- opt$makedb_file
v_germline_file <- opt$v_germline
# Extract the output name from the makedb file (remove extension)
output_name <- tools::file_path_sans_ext(basename(makedb_file))

# Read input data
data <- read.delim(makedb_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
v_germline <- readIgFasta(v_germline_file)

##### Pre filter and collapse identical sequences #####

# Extract V start and sequence
data[["v_start"]] <- stringi::stri_locate_first_regex(data[["sequence_alignment"]], "[ATCG]")
data[["v_seq"]] <- sapply(seq_len(nrow(data)), function(i) substr(data[["sequence_alignment"]][i], 1, data[["v_germline_end"]][i]))

# Calculate mutation count in the V region
data[["v_call_single"]] <- sapply(strsplit(data[["v_call"]], ","), "[[", 1)
data[["v_germline"]] <- sapply(data[["v_call_single"]], function(x) v_germline[[x]])
mutations <- piglet::allele_diff_indices_parallel2(germs = data[["v_germline"]], inputs = data[["v_seq"]])
data[["v_mut"]] <- sapply(seq_len(nrow(data)), function(i) {
  allele <- data[["v_call_single"]][i]
  idx <- mutations[[i]]
  v_min <- min(data[["v_start"]][grep(allele, data[["v_call"]], fixed = TRUE)]) + 5
  sum(idx > v_min & idx <= 316) <= 3
})

# Filter out sequences with mutations in the V region
data <- data[data$v_mut == TRUE, ]

# Collapse identical sequences
if (!"consensus_count" %in% names(data)) {
  data <- data %>%
    mutate(
      consensus_count = ifelse("reads" %in% names(.), reads, 1),
      duplicate_count = ifelse("templates" %in% names(.), templates, 1)
    )
} else {
  data <- data %>%
    mutate(
      consensus_count = ifelse(is.na(consensus_count), 1, consensus_count),
      duplicate_count = ifelse(is.na(duplicate_count), 1, duplicate_count)
    )
}

data <- data %>% select(sequence_id, sequence, sequence_alignment, v_call, d_call, j_call, 
consensus_count, duplicate_count, v_germline_end,
d_sequence_start, d_sequence_end,
d_germline_start, d_germline_end)
data[["sequence_vdj"]] <- gsub("[.]", "", data[["sequence_alignment"]])

data <- data %>% dplyr::group_by(sequence_vdj) %>%
  dplyr::summarise(
    sequence_id = sequence_id[which.max(consensus_count)],
    sequence = sequence[which.max(consensus_count)],
    sequence_alignment = sequence_alignment[which.max(consensus_count)],
    v_call = v_call[which.max(consensus_count)],
    d_call = d_call[which.max(consensus_count)],
    j_call = j_call[which.max(consensus_count)],
    v_germline_end = v_germline_end[which.max(consensus_count)],
    d_sequence_start = d_sequence_start[which.max(consensus_count)],
    d_sequence_end = d_sequence_end[which.max(consensus_count)],
    d_germline_start = d_germline_start[which.max(consensus_count)],
    d_germline_end = d_germline_end[which.max(consensus_count)],
    consensus_count = sum(consensus_count),
    duplicate_count = sum(duplicate_count),
    .groups = "drop")

# Write collapsed data to file
write.table(data, file = paste0(output_name, "_collapsed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Generate FASTA output
seq_names <- apply(data[, c("sequence_id","consensus_count","duplicate_count")], 1, function(x) paste0(names(x), "=", x, collapse = "|"))
seq_names <- gsub('sequence_id=', '', seq_names, fixed = TRUE)
writeFasta(setNames(as.list(data[["sequence"]]), seq_names), file = paste0(output_name, "_collapsed.fasta"))

# End of script
