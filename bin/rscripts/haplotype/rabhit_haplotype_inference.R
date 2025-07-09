#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tigger)
  library(data.table)
  library(rabhit)
  library(alakazam)
})

# ---- Command-line arguments ----
option_list <- list(
  make_option("--input", type="character", help="Input AIRR TSV file"),
  make_option("--v_germline", type="character", help="Path to V germline FASTA"),
  make_option("--d_germline", type="character", help="Path to D germline FASTA"),
  make_option("--prefix", type="character", help="Prefix name for output files", default="rabhit")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Logging function ----
log <- function(msg) {
  cat(paste(Sys.time(), ">", msg), sep = "\n")
}

# ---- Start ----
log("START haplotype inference")
log(paste("Input file:", opt$input))

# Read AIRR data
log("Reading AIRR file...")
data <- fread(opt$input, data.table = FALSE)
log(paste("Sequences loaded:", nrow(data)))

# Read germlines
log("Reading germlines (if provided)...")
v_germline_db <- if (opt$v_germline != "") readIgFasta(opt$v_germline) else NA
d_germline_db <- if (opt$d_germline != "") readIgFasta(opt$d_germline) else NA
prefix <- opt$prefix
# Compute deletion profile
log("Running deletionsByBinom...")
binom_del <- deletionsByBinom(data, chain = "IGH")
outfile_del <- paste0(prefix, "_binomDel.tsv")

write.table(binom_del, file = outfile_del, sep = '\t', row.names = FALSE, quote = TRUE)
log(paste("Deletion report written to:", outfile_del))

# Haplotype inference
log("Checking haplotypable genes...")
outfile_haplotype_prefix <- paste0(prefix, "_gene-")
genes_haplotype <- c('IGHJ6', 'IGHD2-21', 'IGHD2-8')

for (gene in genes_haplotype) {
  CALL <- paste0(tolower(substr(gene, 4, 4)), "_call")

  if (gene == 'IGHJ6') {
    CALL <- 'j_call'
    toHap_GERM <- c(v_germline_db, d_germline_db)
    toHap_col <- c('v_call', 'd_call')
  } else {
    toHap_GERM <- c(v_germline_db)
    toHap_col <- c('v_call')
  }

  allele_fractions <- grep(gene, grep(',', data[[CALL]], invert = TRUE, value = TRUE), value = TRUE)
  bool <- sum(table(allele_fractions) / length(allele_fractions) >= 0.3) == 2 &&
          length(names(table(allele_fractions))) >= 2

  if (bool) {
    names_ <- names(table(allele_fractions)[table(allele_fractions) / length(allele_fractions) >= 0.3])
    alleles <- paste0(sapply(names_, function(x) strsplit(x, '[*]')[[1]][2]), collapse = '_')

    log(paste("Inferring haplotype for", gene, "with alleles:", alleles))

    haplo <- createFullHaplotype(
      data,
      toHap_col = toHap_col,
      hapBy_col = CALL,
      hapBy = gene,
      toHap_GERM = toHap_GERM,
      deleted_genes = binom_del,
      chain = "IGH"
    )

    out_haplo <- paste0(outfile_haplotype_prefix, gene, '-', alleles, "_haplotype.tsv")
    write.table(haplo, file = out_haplo, sep = '\t', row.names = FALSE, quote = TRUE)
    log(paste("Haplotype report written to:", out_haplo))
  }
}

log("END haplotype inference")
