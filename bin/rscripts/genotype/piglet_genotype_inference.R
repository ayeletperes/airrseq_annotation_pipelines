#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(piglet)
  library(tigger)
  library(data.table)
  library(dplyr)
})

# ---- Parse command-line arguments ----
option_list <- list(
  make_option("--input", type = "character"),
  make_option("--germline", type = "character"),
  make_option("--thresholds", type = "character"),
  make_option("--call", type = "character"),
  make_option("--find_unmutated", type = "character"),
  make_option("--single_assignment", type = "character"),
  make_option("--z_threshold", type = "numeric")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Logging ----
log <- function(msg) cat(paste(Sys.time(), ">", msg), sep = "\n")
log("START Piglet genotype inference")

# ---- Load inputs ----
data <- fread(opt$input, data.table = FALSE)
log(paste("Input rows:", nrow(data)))

if (opt$thresholds == "") {
  stop("No allele threshold table supplied")
} else {
  allele_threshold_table <- fread(opt$thresholds)
}

germline_db <- if (opt$germline != "") readIgFasta(opt$germline) else NA

# ---- Expand threshold table for novel alleles ----
if (length(germline_db) > 0 && any(grepl("_", names(germline_db)))) {
  alleles <- grep("_", names(germline_db), value = TRUE)
  for (a in alleles) {
    a_split <- unlist(strsplit(a, "_"))
    base_allele <- a_split[1]
    snps <- paste0(a_split[-1], collapse = "_")
    base_threshold <- allele_threshold_table[asc_allele == base_allele, ]
    if (nrow(base_threshold) > 0) {
      base_threshold$asc_allele <- a
      base_threshold$allele <- paste0(base_threshold$allele, "_", snps)
      allele_threshold_table <- rbind(allele_threshold_table, base_threshold)
    }
  }
}

# ---- Run Piglet inference ----
find_unmutated <- opt$find_unmutated == "true"
single_assignment <- opt$single_assignment == "true"
call_col <- opt$call
z_thresh <- opt$z_threshold

log("Running Piglet::inferGenotypeAllele")
geno <- piglet::inferGenotypeAllele(
  data,
  allele_threshold_table = allele_threshold_table,
  germline_db = germline_db,
  find_unmutated = find_unmutated,
  call = call_col,
  single_assignment = single_assignment
)

## --- Filter by z-score for personal reference --- 
geno_alleles <- geno[geno$z_score >= z_thresh, ]

# ---- Build personal reference ----
NOTGENO.IND <- !(sapply(strsplit(names(germline_db), "*", fixed = TRUE), `[`, 1) %in% geno$gene)
germline_db_new <- germline_db[NOTGENO.IND]

for (i in seq_len(nrow(geno))) {
  allele <- geno$allele[i]
  IND <- names(germline_db) %in% allele
  germline_db_new <- c(germline_db_new, germline_db[IND])
}

ref_out <- paste0(call_col, "_personal_reference.fasta")
writeFasta(germline_db_new, file = ref_out)
log(paste("Personal reference written to:", ref_out))

# ---- Transform genotype output to match VDJbase format
geno[, allele_post_star := sub(".*[*]", "", allele)]
geno[, genotyped_allele_temp := ifelse(z_score > 0, allele_post_star, NA)]

# Summarize per gene
geno <- geno[, .(
  gene = gene[1],
  alleles = paste(allele_post_star, collapse=","),
  sample = paste(unique(sample), collapse=","),
  threshold = paste(threshold, collapse=","),
  counts = paste(count, collapse=","),
  total = sum(count),
  note = "",
  z_score = paste(z_score, collapse=","),
  genotyped_alleles = paste(na.omit(genotyped_allele_temp), collapse=",")
), by = gene]

# Write full report
report_path <- paste0(call_col, "_genotype_report.tsv")
write.table(geno, file = report_path, row.names = FALSE, sep = "\t")
log(paste("Genotype report written to:", report_path))

log("END Piglet genotype inference")
