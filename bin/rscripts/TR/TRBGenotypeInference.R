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
  make_option(c("-v", "--v_germline"), type="character", help="Path to the V germline FASTA file", metavar="FILE"),
  make_option(c("-d", "--d_germline"), type="character", help="Path to the D germline FASTA file", metavar="FILE"),
  make_option(c("-j", "--j_germline"), type="character", help="Path to the J germline FASTA file", metavar="FILE"),
  make_option(c("-s", "--max_snp_position"), type="numeric", help="The max snp position for counting mutation to filter", default=316),
  make_option(c("-c", "--min_consensus_count"), type="numeric", help="The minimal consensus count to filter", default=1)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$makedb_file) || is.null(opt$v_germline) || is.null(opt$d_germline) || is.null(opt$j_germline)) {
  print_help(opt_parser)
  quit(status=1)
}

# Assign variables
makedb_file <- opt$makedb_file
v_germline_file <- opt$v_germline
d_germline_file <- opt$d_germline
j_germline_file <- opt$j_germline
max_snp_position <- opt$max_snp_position
min_consensus_count <- opt$min_consensus_count

# Read input data
data <- read.delim(makedb_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
v_germline <- readIgFasta(v_germline_file)
d_germline <- readIgFasta(d_germline_file)
j_germline <- readIgFasta(j_germline_file)

## genotype inference
# Filter by consensus count
data <- data[data[["consensus_count"]] >= min_consensus_count, ]

# Calculate mutation count in the V region and filter to zero mutations
data[["v_seq"]] <- substr(data[["sequence_alignment"]], 1, sapply(as.numeric(data[["v_germline_end"]]), min, max_snp_position)) 
data[["v_mut"]] <- sapply(tigger::getMutCount(data[["v_seq"]], data[["v_call"]], germline_db = v_germline), function(x){x[[1]]}) 
data <- data[data[["v_mut"]] <= 1, ]


DATA_V_SA <- data[!grepl(pattern = ',', data[["v_call"]]), ]

geno_BV <- inferGenotypeBayesian(DATA_V_SA, germline_db = v_germline, find_unmutated = FALSE, novel = new_novel_df_H, v_call = 'v_call')
names(geno_BV) <- names(geno_BV)
geno_BV[["genotyped_alleles"]] <- apply(geno_BV[, c(2, 6:9)], 1, function(y){m <- which.max(as.numeric(y[2:5]));paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")})

DATA_D_geno <- data[(!grepl(pattern = ',', data[["d_call"]]) & data[["d_call"]] != 'None') & (data[["d_sequence_end"]] - data[["d_sequence_start"]] >= 8), ]
DATA_D_geno <- DATA_D_geno[complete.cases(DATA_D_geno[["sequence_id"]]), ]

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

# filter by zero mutations over the D segment
DATA_D_geno[["mut_d"]] <- unlist(lapply(seq_len(nrow(DATA_D_geno)), function(i) {
  mut <- 0
  row_seq <- unlist(strsplit(DATA_D_geno[["d_seq"]][[i]], ""))
  allele_seq <- unlist(strsplit(d_germline[[DATA_D_geno[["d_call"]][[i]]]], ""))
  for (pos in DATA_D_geno[["d_germline_start"]][[i]]:DATA_D_geno[["d_germline_end"]][[i]]) {
    if (row_seq[pos - (DATA_D_geno[["d_germline_start"]][[i]] - 1)] != allele_seq[pos]) {
      mut <- mut + 1
    }
  }
  mut
}))

DATA_D_geno <- DATA_D_geno[DATA_D_geno[["mut_d"]] == 0, ]

geno_BD <- inferGenotypeBayesian(DATA_D_geno, find_unmutated = FALSE, germline_db = d_germline, v_call = 'd_call')
geno_BD[["genotyped_alleles"]] <- apply(geno_BD[, c(2, 6:9)], 1, function(y){m <- which.max(as.numeric(y[2:3]));paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")})

D2_total <- nrow(DATA_D_geno[grepl("TRBD2", DATA_D_geno[["d_call"]]), ])
D2_01_count <- nrow(DATA_D_geno[DATA_D_geno[["d_call"]] == "TRBD2*01", ])
D2_01_freq <- D2_01_count / D2_total

if (D2_01_freq < 0.2066) {
  geno_BD[["genotyped_alleles"]][geno_BD[["gene"]] == "TRBD2"] <- "02"
} else if (D2_01_freq > 0.8969) {
  geno_BD[["genotyped_alleles"]][geno_BD[["gene"]] == "TRBD2"] <- "01"
} else if (geno_BD[["genotyped_alleles"]][geno_BD[["gene"]] == "TRBD2"] == "01") {
  geno_BD[["genotyped_alleles"]][geno_BD[["gene"]] == "TRBD2"] <- "01,02"
}

DATA_J_SA <- data[!grepl(pattern = ',', data[["j_call"]]), ]
geno_BJ <- inferGenotypeBayesian(data, germline_db = j_germline, find_unmutated = FALSE, v_call = 'j_call')
geno_BJ[["genotyped_alleles"]] <- apply(geno_BJ[, c(2, 6:9)], 1, function(y){m <- which.max(as.numeric(y[2:3]));paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")})


## Remove from v_germline irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(v_germline),'*',fixed=T),'[',1) %in%  geno_BV[["gene"]])
TRBV_GERM.NEW <- v_germline[NOTGENO.IND]

for(i in seq_len(nrow(geno_BV))){
  gene <- geno_BV[["gene"]][i]
  
  alleles <- geno_BV[["genotyped_alleles"]][i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(v_germline) %in%  paste(gene,alleles,sep='*')
  TRBV_GERM.NEW <- c(TRBV_GERM.NEW,v_germline[IND])
}


## Remove from d_germline irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(d_germline),'*',fixed=T),'[',1) %in%  geno_BD[["gene"]])
TRBD_GERM.NEW <- d_germline[NOTGENO.IND]

for(i in seq_len(nrow(geno_BD))){
  gene <- geno_BD[["gene"]][i]
  alleles <- geno_BD[["genotyped_alleles"]][i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(d_germline) %in%  paste(gene,alleles,sep='*')
  TRBD_GERM.NEW <- c(TRBD_GERM.NEW,d_germline[IND])
}

## Remove from j_germline irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(j_germline),'*',fixed=T),'[',1) %in%  geno_BJ[["gene"]])
TRBJ_GERM.NEW <- j_germline[NOTGENO.IND]

for(i in seq_len(nrow(geno_BJ))){
  gene <- geno_BJ[["gene"]][i]
  alleles <- geno_BJ[["genotyped_alleles"]][i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(j_germline) %in%  paste(gene,alleles,sep='*')
  TRBJ_GERM.NEW <- c(TRBJ_GERM.NEW,j_germline[IND])
}


### CHECK IF THE REPLACEMENT IS CORRECT

## Combine the genotyped and others and write to a fasta file for reference

writeFasta(TRBV_GERM.NEW, file = "v_personal.fasta")
writeFasta(TRBD_GERM.NEW, file = "d_personal.fasta")
writeFasta(TRBJ_GERM.NEW, file = "j_personal.fasta")

## save the genotype data
write.table(geno_BV, file = "v_genotype.tsv", quote = F, row.names = F, sep = "\t")
write.table(geno_BD, file = "d_genotype.tsv", quote = F, row.names = F, sep = "\t")
write.table(geno_BJ, file = "j_genotype.tsv", quote = F, row.names = F, sep = "\t")
# End of script
