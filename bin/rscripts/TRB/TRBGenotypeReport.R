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
  make_option(c("-v", "--v_genotype"), type="character", help="Path to the V genotype file", metavar="FILE"),
  make_option(c("-d", "--d_genotype"), type="character", help="Path to the D genotype file", metavar="FILE"),
  make_option(c("-j", "--j_genotype"), type="character", help="Path to the J genotype file", metavar="FILE"),
  make_option(c("-o", "--outname"), type="character", help="Output file name", default="TRB")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$makedb_file) || is.null(opt$v_genotype) || is.null(opt$d_genotype) || is.null(opt$j_genotype)) {
  print_help(opt_parser)
  quit(status=1)
}

# Assign variables
makedb_file <- opt$makedb_file
v_genotype_file <- opt$v_genotype
d_genotype_file <- opt$d_genotype
j_genotype_file <- opt$j_genotype
outname <- opt$outname
# Read input data
data <- read.delim(makedb_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
geno_BV <- read.delim(v_genotype_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
geno_BD <- read.delim(d_genotype_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
geno_BJ <- read.delim(j_genotype_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

v_table <- data[!grepl(",", data[["v_call"]]), ]
v_table <- v_table %>% dplyr::group_by(v_call) %>% dplyr::summarise(N = n())
v_table[["gene"]] <- sapply(strsplit(v_table[["v_call"]], "*", fixed = TRUE), "[", 1)
v_table[["allele"]] <- sapply(strsplit(v_table[["v_call"]], "*", fixed = TRUE), "[", 2)
geno_BV[["Freq_by_Seq"]] <- unlist(lapply(geno_BV[["gene"]], function(g) {
  if (!g %in% v_table[["gene"]]) {
    return("0")
  }
  counts <- v_table[["N"]][v_table[["gene"]] == g]
  counts <- sort(counts, decreasing = TRUE)
  return(paste(counts, collapse = ";"))
}))

d_table <- data[(!grepl(pattern = ',', data[["d_call"]]) & data[["d_call"]] != 'None') & (data[["d_sequence_end"]] - data[["d_sequence_start"]] >= 8), ]
d_table <- d_table[complete.cases(d_table[["sequence_id"]]), ]
d_table <- d_table %>% dplyr::group_by(d_call) %>% dplyr::summarise(N = n())
d_table[["gene"]] <- sapply(strsplit(d_table[["d_call"]], "*", fixed = TRUE), "[", 1)
d_table[["allele"]] <- sapply(strsplit(d_table[["d_call"]], "*", fixed = TRUE), "[", 2)
geno_BD[["Freq_by_Seq"]] <- unlist(lapply(geno_BD[["gene"]], function(g) {
  if (!g %in% d_table[["gene"]]) {
    return("0")
  }
  counts <- d_table[["N"]][d_table[["gene"]] == g]
  counts <- sort(counts, decreasing = TRUE)
  return(paste(counts, collapse = ";"))
}))

j_table <- data[!grepl(",", data[["j_call"]]), ]
j_table <- j_table %>% dplyr::group_by(j_call) %>% dplyr::summarise(N = n())
j_table[["gene"]] <- sapply(strsplit(j_table[["j_call"]], "*", fixed = TRUE), "[", 1)
j_table[["allele"]] <- sapply(strsplit(j_table[["j_call"]], "*", fixed = TRUE), "[", 2)
geno_BJ[["Freq_by_Seq"]] <- unlist(lapply(geno_BJ[["gene"]], function(g) {
  if (!g %in% j_table[["gene"]]) {
    return("0")
  }
  counts <- j_table[["N"]][j_table[["gene"]] == g]
  counts <- sort(counts, decreasing = TRUE)
  return(paste(counts, collapse = ";"))
}))

geno <- rbind(geno_BV, geno_BD, geno_BJ)
geno[["Freq_by_Clone"]] <- geno[["Freq_by_Seq"]]
write.table(geno, file = paste0(outname,"_genotype.tsv"), quote = FALSE, row.names = F, sep = "\t")