#!/usr/bin/env Rscript

# Load required libraries
suppressMessages(library(piglet))
suppressMessages(library(tigger))
suppressMessages(library(dplyr))
library(stringr)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("-i", "--i"), type="character", help="Prefix name for output files", default=""),
  make_option(c("-m", "--makedb_file"), type="character", help="Path to the makedb file (input TSV)", metavar="FILE"),
  make_option(c("-a", "--v_germline"), type="character", help="Path to the V germline FASTA file", metavar="FILE"),
  make_option(c("-b", "--v_genotype"), type="character", help="Path to the V genotype TSV file", metavar="FILE"),
  make_option(c("-u", "--gene_usages_file"), type="numeric", help="Gene usage file, for deletion inference.", metavar="FILE"),
  make_option(c("-c", "--min_consensus_count"), type="numeric", help="The minimal consensus count to filter", default=1)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$makedb_file) || is.null(opt$v_germline) || is.null(opt$v_genotype) || is.null(opt$gene_usages_file)) {
  print_help(opt_parser)
  quit(status=1)
}

# Assign variables
makedb_file <- opt$makedb_file
v_germline_file <- opt$v_germline
v_genotype_file <- opt$v_genotype
gene_usages_file <- opt$gene_usages_file
min_consensus_count <- opt$min_consensus_count
prefix <- opt$i

# Read input data
DATA <- read.delim(makedb_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
TRBV_GERM <- readIgFasta(v_germline_file)
geno_BV <- read.delim(v_genotype_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_usages <- read.csv(gene_usages_file, sep = "\t")

if ("consensus_count" %in% colnames(DATA)) {
  DATA <- DATA[DATA[["consensus_count"]] >= min_consensus_count, ]
} else {
  print("consensus_count column not found in DATA")
}

cutoff <- 0.0005
p_val_cutoff <- 0.001

DATA[["v_gene"]] <- unlist(lapply(DATA[["v_call"]], function(Vcall){
  assignments <- unlist(strsplit(Vcall, ",", fixed = TRUE))
  v_genes <- unique(sapply(strsplit(assignments, "*", fixed = TRUE), "[", 1))
  v_genes <- v_genes[order(v_genes)]
  paste(v_genes, collapse = ",")
}))

DATA <- DATA[!grepl(",", DATA[["v_gene"]]), ]
DATA <- DATA[grepl("TRBV", DATA[["v_gene"]]), ]

gene_usages[["N"]] <- unlist(lapply(gene_usages[["GENE"]], function(gene){nrow(DATA[DATA[["v_gene"]] == gene, ])}))
gene_usages[["TOTAL"]] <- nrow(DATA)
gene_usages[["USAGE"]] <- gene_usages[["N"]] / gene_usages[["TOTAL"]]

gene_usages <- gene_usages[gene_usages[["MIN_FREQ"]] != Inf & gene_usages[["AVG_USAGE"]] > 1.5 * cutoff, ]

# calculate the p value for each gene in case that the gene is deleted by binom test
gene_usages[["PVAL"]] <- sapply(1:nrow(gene_usages), function(i) {
  if (gene_usages[["USAGE"]][i] < cutoff) {
    return(binom.test(x = gene_usages[["N"]][i], n = gene_usages[["TOTAL"]][i], p = gene_usages[["MIN_FREQ"]][i])$p.value)
  } else {
    return(1)
  }
})

# Detect according to the p values if there are deleted genes
gene_usages[["DELETED"]] <- sapply(1:nrow(gene_usages), function(i) {
  if (gene_usages[["PVAL"]][i] <= p_val_cutoff) {
    if ((gene_usages[["USAGE"]][i] < cutoff) & gene_usages[["MIN_FREQ"]][i] != Inf) {
      return(TRUE)
    }
  }
  return(FALSE)
})

gene_usages <- gene_usages[gene_usages[["DELETED"]], ]

if (nrow(gene_usages) > 0) {
  deleted_genes <- gene_usages[["GENE"]]
  geno_BV <- geno_BV[!geno_BV[["gene"]] %in% deleted_genes, ]
  
  for (gene in deleted_genes) {
    geno_BV[nrow(geno_BV) + 1, ] <- c(gene, NA, NA, NA, NA, NA, NA, NA, NA, 1000, "Deletion")
  }
}

## save the deletion data
write.table(geno_BV, file = paste0(prefix,"_v_genotype_deletion.tsv"), quote = F, row.names = F, sep = "\t")
# End of script
