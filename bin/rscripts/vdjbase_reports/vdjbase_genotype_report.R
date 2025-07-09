#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(alakazam)
  library(plyr)
})

# ---- Command-line arguments ----
option_list <- list(
  make_option("--initial_reads", type="character"),
  make_option("--personal_reads", type="character"),
  make_option("--v_genotype", type="character"),
  make_option("--d_genotype", type="character"),
  make_option("--j_genotype", type="character"),
  make_option("--outname", type="character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Logging function ----
log <- function(msg) {
  cat(paste(Sys.time(), ">", msg), sep = "\n")
}

log("START VDJbase genotype report")

getFreq <- function(data, call = "v_call") {
  table(grep(",", data[[call]][data[[call]] != ""], invert = TRUE, value = TRUE))
}

addFreqInfo <- function(tab, gene, alleles) {
  paste0(tab[paste0(gene, "*", unlist(strsplit(alleles, ',')))], collapse = ";")
}

# ---- Load data ----
log("Reading input tables...")
data_initial <- read_tsv(opt$initial_reads, col_select = c("sequence_id", "v_call", "d_call", "j_call"))
data_genotyped <- read_tsv(opt$personal_reads, col_select = c("sequence_id", "v_call", "d_call", "j_call"))

# ---- Match sequence IDs ----
data_initial <- data_initial[data_initial$sequence_id %in% data_genotyped$sequence_id, ]
data_genotyped <- data_genotyped[data_genotyped$sequence_id %in% data_initial$sequence_id, ]
data_initial <- data_initial[order(data_initial$sequence_id), ]
data_genotyped <- data_genotyped[order(data_genotyped$sequence_id), ]

# ---- Fix mismatches in V calls ----
non_match_v <- which(data_initial$v_call != data_genotyped$v_call)
data_initial$v_call[non_match_v] <- data_genotyped$v_call[non_match_v]

# ---- V genotype frequencies ----
log("Computing v_call frequencies...")
tab_freq_v <- getFreq(data_genotyped, call = "v_call")
tab_clone_v <- getFreq(data_initial, call = "v_call")
tab_clone_v <- tab_clone_v[names(tab_freq_v)]
genoV <- read_tsv(opt$v_genotype, col_types = cols(.default = "c"))
log("Adding v_call frequencies...")
genoV <- genoV %>% dplyr::group_by(gene) %>% dplyr::mutate(
  Freq_by_Clone = addFreqInfo(tab_clone_v, gene, genotyped_alleles),
  Freq_by_Seq = addFreqInfo(tab_freq_v, gene, genotyped_alleles)
)

# ---- J genotype frequencies ----
log("Computing j_call frequencies...")
tab_freq_j <- getFreq(data_genotyped, call = "j_call")
tab_clone_j <- getFreq(data_initial, call = "j_call")
tab_clone_j <- tab_clone_j[names(tab_freq_j)]
genoJ <- read_tsv(opt$j_genotype, col_types = cols(.default = "c"))
genoJ <- genoJ %>% dplyr::group_by(gene) %>% dplyr::mutate(
  Freq_by_Clone = addFreqInfo(tab_clone_j, gene, genotyped_alleles),
  Freq_by_Seq = addFreqInfo(tab_freq_j, gene, genotyped_alleles)
)

# ---- Optional D genotype frequencies ----
if (endsWith(opt$d_genotype, ".tsv")) {
  log("Computing d_call frequencies...")
  tab_freq_d <- getFreq(data_genotyped, call = "d_call")
  tab_clone_d <- getFreq(data_initial, call = "d_call")
  tab_clone_d <- tab_clone_d[names(tab_freq_d)]
  genoD <- read_tsv(opt$d_genotype, col_types = cols(.default = "c"))
  genoD <- genoD %>% dplyr::group_by(gene) %>% dplyr::mutate(
    Freq_by_Clone = addFreqInfo(tab_clone_d, gene, genotyped_alleles),
    Freq_by_Seq = addFreqInfo(tab_freq_d, gene, genotyped_alleles)
  )
  genos <- plyr::rbind.fill(genoV, genoD, genoJ)
} else {
  genos <- plyr::rbind.fill(genoV, genoJ)
}

# ---- Cleanup and write ----
genos$Freq_by_Clone <- gsub("NA", "0", genos$Freq_by_Clone)
genos$Freq_by_Seq <- gsub("NA", "0", genos$Freq_by_Seq)
colnames(genos)[colnames(genos) == "genotyped_alleles"] <- "GENOTYPED_ALLELES"

output_file <- paste0(opt$outname, "_Final_genotype.tsv")
log(paste("Writing final report to:", output_file))
write.table(genos, file = output_file, row.names = FALSE, sep = "\t")

log("END VDJbase genotype report")
