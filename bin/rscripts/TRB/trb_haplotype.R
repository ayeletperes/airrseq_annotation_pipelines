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
  make_option("--genos", type="character", help="Path to genotype file"),
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
#data <- fread(opt$input, data.table = FALSE)
#log(paste("Sequences loaded:", nrow(data)))

# Read germlines
log("Reading germlines (if provided)...")
TRBV_GERM <- readIgFasta(opt$v_germline)
TRBD_GERM <- readIgFasta(opt$d_germline)

geno <- read.delim(opt$genos, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DATA <- read.delim(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
prefix <- opt$prefix

DATA[["subject"]] <- prefix

#hetero_j16 <- "TRBJ1-6" %in% geno[["gene"]]
#if (hetero_j16) {
#  hetero_j16 <- grepl(",", geno[["GENOTYPED_ALLELES"]][geno[["gene"]] == "TRBJ1-6"])
#}

#hetero_d2 <- "TRBD2" %in% geno[["gene"]]
#if (hetero_d2) {
#  hetero_d2 <- grepl(",", geno[["GENOTYPED_ALLELES"]][geno[["gene"]] == "TRBD2"])
#}



hetero_j16 <- "TRBJ1-6" %in% geno[["gene"]] &&
              any(grepl(",", geno[["genotyped_alleles"]][geno[["gene"]] == "TRBJ1-6"]))

hetero_d2 <- "TRBD2" %in% geno[["gene"]] &&
              any(grepl(",", geno[["genotyped_alleles"]][geno[["gene"]] == "TRBD2"]))

# Filtering
DATA <- DATA[!grepl(",", DATA[["v_call"]]), ] # V single assignment

# zero mutations over the V
v_seqs <- sapply(1:nrow(DATA), function(x) substr(DATA[["sequence_alignment"]][x], 1, DATA[["v_germline_end"]][x]))
DATA[["v_mut"]] <- unlist(tigger::getMutCount(v_seqs, DATA[["v_call"]], germline_db = TRBV_GERM))
DATA <- DATA[DATA[["v_mut"]] <= 1, ]

del_genes <- geno[["gene"]][grepl("Del", geno[["genotyped_alleles"]])]

if (hetero_j16) {
  # Only rearrangements with single assignment of TRBJ1-6
  DATA_J <- DATA[grepl("J1-6", DATA[["j_call"]]), ]
  DATA_J <- DATA_J[!grepl(",", DATA_J[["j_call"]]), ]
  DATA_J[["J_SEQ_LENGTH"]] <- DATA_J[["j_sequence_end"]] - DATA_J[["j_sequence_start"]] + 1
  DATA_J <- DATA_J[DATA_J[["J_SEQ_LENGTH"]] > 10, ]
  
  TRBJ1_6_01_REA <- nrow(DATA_J[DATA_J[["j_call"]] == "TRBJ1-6*01", ])
  total <- nrow(DATA_J)
  
  if ((TRBJ1_6_01_REA / total > 0.3) & (TRBJ1_6_01_REA / total < 0.7)) {
    haplo_j1_6 <- createFullHaplotype(DATA_J, toHap_col = "v_call", hapBy_col = "j_call", chain = "TRB", hapBy = "TRBJ1-6", toHap_GERM = TRBV_GERM, rmPseudo = FALSE)
    haplo_j1_6[["TOTAL"]] <- total
    
    if (length(del_genes) > 0) {
      for (gene in del_genes) {
        haplo_j1_6[haplo_j1_6[["gene"]] == gene, c(3:5, 9)] <- c("Del", "Del", "Del", 1000)
      }
    }
    
    write.table(haplo_j1_6, file = paste0(prefix,"_gene-TRBJ1_6_haplotype.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}

if (hetero_d2) {
  # Only rearrangements with single assignment of TRBD2
  DATA[["d_germline_length"]] <- as.integer(DATA[["d_germline_end"]]) - as.integer(DATA[["d_germline_start"]]) + 1
  DATA_D_geno <- DATA[(!grepl(pattern = ',', DATA[["d_call"]]) & DATA[["d_call"]] != 'None') & (DATA[["d_germline_length"]] >= 9), ]
  DATA_D_geno <- DATA_D_geno[complete.cases(DATA_D_geno[["sequence_id"]]), ]
  DATA_D_geno <- DATA_D_geno[grepl("D2", DATA_D_geno[["d_call"]]), ]
  
  # filter by zero mutations over the D segment
  # extract d sequence in the direct orientation
  DATA_D_reg <- DATA_D_geno[DATA_D_geno[["d_germline_start"]] < DATA_D_geno[["d_germline_end"]], ]
  DATA_D_reg[["d_seq"]] <- substr(DATA_D_reg[["sequence"]], DATA_D_reg[["d_sequence_start"]], DATA_D_reg[["d_sequence_end"]])
  
  # extract convert d sequence in the inverted orientation to the direct orientation
  DATA_D_inv <- DATA_D_geno[DATA_D_geno[["d_germline_start"]] > DATA_D_geno[["d_germline_end"]], ]
  DATA_D_inv[["d_seq"]] <- substr(DATA_D_inv[["sequence"]], DATA_D_inv[["d_sequence_start"]], DATA_D_inv[["d_sequence_end"]])
  DATA_D_inv[["d_seq"]] <- stringi::stri_reverse(DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- gsub("A", "t", DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- gsub("T", "a", DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- gsub("G", "c", DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- gsub("C", "g", DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- toupper(DATA_D_inv[["d_seq"]])
  
  d_germ_end <- DATA_D_inv[["d_germline_start"]]
  DATA_D_inv[["d_germline_start"]] <- DATA_D_inv[["d_germline_end"]]
  DATA_D_inv[["d_germline_end"]] <- d_germ_end
  
  DATA_D_geno <- rbind(DATA_D_reg, DATA_D_inv)
  
  DATA_D_geno[["mut_d"]] <- unlist(lapply(1:nrow(DATA_D_geno), function(i) {
    mut <- 0
    row_seq <- unlist(strsplit(DATA_D_geno[["d_seq"]][[i]], ""))
    allele_seq <- unlist(strsplit(TRBD_GERM[[DATA_D_geno[["d_call"]][[i]]]], ""))
    for (pos in DATA_D_geno[["d_germline_start"]][[i]]:(DATA_D_geno[["d_germline_end"]][[i]])) {
      if (row_seq[pos - (DATA_D_geno[["d_germline_start"]][[i]] - 1)] != allele_seq[pos]) {
        mut <- mut + 1
      }
    }
    mut
  }))
  
  DATA_D_geno <- DATA_D_geno[DATA_D_geno[["mut_d"]] == 0, ]
  TRBD2_01_REA <- nrow(DATA_D_geno[DATA_D_geno[["d_call"]] == "TRBD2*01", ])
  total <- nrow(DATA_D_geno)
  
  if ((TRBD2_01_REA / total > 0.3) & (TRBD2_01_REA / total < 0.7)) {
    haplo_d2 <- createFullHaplotype(DATA_D_geno, toHap_col = "v_call", hapBy_col = "d_call", chain = "TRB", hapBy = "TRBD2", toHap_GERM = TRBV_GERM, rmPseudo = FALSE, kThreshDel = 5)
    haplo_d2[["TOTAL"]] <- total
    
    if (length(del_genes) > 0) {
      for (gene in del_genes) {
        haplo_d2[haplo_d2[["gene"]] == gene, c(3:5, 9)] <- c("Del", "Del", "Del", 1000)
      }
    }
    
    write.table(haplo_d2, file =  paste0(prefix,"_gene-TRBD2_1_2_haplotype.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}


log("END haplotype inference")
