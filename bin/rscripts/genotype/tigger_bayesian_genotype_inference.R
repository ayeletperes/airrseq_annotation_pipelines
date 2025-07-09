#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tigger)
  library(data.table)
  library(stringr)
})

# ---- Parse command-line arguments ----
option_list <- list(
  make_option("--input", type="character", help="Input AIRR TSV file"),
  make_option("--germline", type="character", help="Path to germline FASTA file"),
  make_option("--call", type="character", help="Column name of segment call (e.g., v_call)"),
  make_option("--find_unmutated", type="character", help="true/false"),
  make_option("--single_assignments", type="character", help="true/false")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Prepare log file

log <- function(msg) {
  cat(paste(Sys.time(), ">", msg), sep = "\n")
}

# ---- Start logging
log("START genotype inference")
log(paste("Method: TIgGER Bayesian"))
log(paste("Input file:", opt$input))
log(paste("Call column:", opt$call))

# ---- Read input data
data <- fread(opt$input, data.table = FALSE)
log(paste("Input sequences:", nrow(data)))

find_unmutated <- opt$find_unmutated == "true"
single_assignments <- opt$single_assignments == "true"
call_col <- opt$call

if (single_assignments) {
  before_filter <- nrow(data)
  data <- data[!grepl(",", data[[call_col]]), ]
  log(paste("Removed", before_filter - nrow(data), "sequences with multiple assignments"))
}

data <- data[!is.na(data[[call_col]]), ]
log(paste("Remaining sequences after NA filtering:", nrow(data)))

germline_db <- if (opt$germline != "") readIgFasta(opt$germline) else NA

# ---- TiGGer priors
priors <- list(
  v_call = c(0.6, 0.4, 0.4, 0.35, 0.25, 0.25, 0.25, 0.25, 0.25),
  d_call = c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0),
  j_call = c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0)
)

log("Running inferGenotypeBayesian...")

geno <- inferGenotypeBayesian(
    data,
    find_unmutated = find_unmutated,
    germline_db = germline_db,
    v_call = call_col,
    priors = priors[[call_col]]
)


GENOTYPED_ALLELES <- function(y) {
m <- which.max(as.numeric(y[2:5]))
paste0(unlist(strsplit(y[1], ","))[1:m], collapse = ",")
}

geno$genotyped_alleles <- apply(geno[, c(2, 6:9)], 1, GENOTYPED_ALLELES)

out_report <- paste0(call_col, "_genotype_report.tsv")
write.table(geno, file = out_report, sep = "\t", quote = FALSE, row.names = FALSE)
log(paste("Genotype report written to:", out_report))

# ---- Build personal reference
NOTGENO.IND <- !(sapply(strsplit(names(germline_db), "*", fixed = TRUE), `[`, 1) %in% geno$gene)
germline_db_new <- germline_db[NOTGENO.IND]

for (i in seq_len(nrow(geno))) {
gene <- geno[i, "gene"]
alleles <- if (geno[i, "genotyped_alleles"] == "") geno[i, "alleles"] else geno[i, "genotyped_alleles"]
alleles <- unlist(strsplit(alleles, ","))
IND <- names(germline_db) %in% paste(gene, alleles, sep = "*")
germline_db_new <- c(germline_db_new, germline_db[IND])
}

ref_out <- paste0(call_col, "_personal_reference.fasta")
writeFasta(germline_db_new, file = ref_out)

log(paste("Personal reference written to:", ref_out))

log("END genotype inference")

