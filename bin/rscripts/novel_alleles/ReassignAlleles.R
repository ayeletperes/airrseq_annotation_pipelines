#!/usr/bin/env Rscript
suppressMessages({
    library(optparse)
    require(tigger)
})

option_list <- list(
  make_option("--input", type = "character", dest = "input"),
  make_option("--out_name", type = "character", dest = "out_name"),
  make_option("--germline_file", type = "character", dest = "germline_file"),
  make_option("--novel_file", type = "character", dest = "novel_file"),
  make_option("--call_column", type = "character", default="v_call", dest = "call_column")
)

opt <- parse_args(OptionParser(option_list = option_list))

log <- function(msg) {
  cat(paste(Sys.time(), ">", msg), sep = "\n")
}

log("Starting Undocumented Allele Reassignment")

input_file <- opt$input
out_prefix <- opt$out_name
germline_file <- opt$germline_file
novel_file <- opt$novel_file
call_column <- opt$call_column

log("Reading input data and Undocumented Allele")

data <- data.table::fread(input_file, stringsAsFactors = FALSE, data.table = FALSE)
novel <- data.table::fread(novel_file, stringsAsFactors = FALSE, data.table = FALSE)
germline <- readIgFasta(germline_file)

log("Getting novel alleles")
novel_alleles <- novel[['polymorphism_call']]
novel_genes <- novel[['gene']]
## get all the assignments that has the novel genes
novel_genes_regex <- paste0(novel_genes, collapse = "|")
novel_assignments_idx <- grep(novel_genes_regex, data[[call_column]])

log("Extracting potential sequences")
novel_assignments <- data[[call_column]][novel_assignments_idx]
novel_assignments <- unique(unlist(strsplit(unique(novel_assignments), ",")))
germline_novel <- germline[c(novel_alleles, novel_assignments)]

data_sub <- data[novel_assignments_idx,c(call_column, "sequence_alignment")]

log("Running reassignAlleles")
## run tigger reassignAlleles
reassigned_db <- reassignAlleles(data_sub, germline_novel, v_call=call_column, seq="sequence_alignment")

## count the number of reassignmetns
num_reassignments <- sum(reassigned_db[['v_call_genotyped']] != reassigned_db[[call_column]], na.rm = TRUE)
log(paste("Number of reassignments:", num_reassignments))

## replace the input column column with the new 'v_call_genotyped' column  and remove `v_call_genotyped`
data[[call_column]][novel_assignments_idx] <- reassigned_db[['v_call_genotyped']]

log("Writing output")
data.table::fwrite(data, file = sprintf('%s_%s_reassigned-pass.tsv',out_prefix, call_column), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

