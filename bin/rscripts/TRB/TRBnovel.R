#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(optparse)
})

# ---- Command-line arguments ----
#option_list <- list(
#  make_option("--input", type="character", help="Input AIRR TSV file"),
#  make_option("-v", "--v_germline", type="character", help="Path to the V germline FASTA file", metavar="FILE")
#)


option_list <- list(
  make_option(c("-m", "--input"), type="character", help="Path to the makedb file (input TSV)", metavar="FILE"),
  make_option(c("-v", "--v_germline"), type="character", help="Path to the V germline FASTA file", metavar="FILE"))
  
opt <- parse_args(OptionParser(option_list = option_list))

data <- read.delim(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
TRBV_GERM <- tigger::readIgFasta(opt$v_germline)


DATA <- data

position <- 316

#DATA <- DATA[grepl("^[A-Z]",DATA$sequence_alignment), ]

#DATA <- DATA[!grepl("^\\.{10,}", DATA$sequence_alignment), ]

n_before <- nrow(DATA)
# Filter out rows starting with at least 10 dots
DATA <- DATA[!grepl("^\\.{10,}", DATA$sequence_alignment), ]
# Count rows after filtering
n_after <- nrow(DATA)
# Calculate and print fraction deleted
fraction_deleted <- (n_before - n_after) / n_before
cat("Fraction of rows deleted:", fraction_deleted, "\n")


DATA[["sequence_v1"]] <- mapply(function(seq, start, end) {
  #substr(seq, start, min(end, position))
  substr(seq, start, position)
  #s <- substr(seq, start, pmin(end, position))
  #gsub("[.]", "", s)
}, DATA[["sequence_alignment"]], DATA[["v_germline_start"]], DATA[["v_germline_end"]])


DATA[["germline_v"]] <- mapply(function(seq, start, end) {
  substr(seq, start, end)
  #s <- substr(seq, start, pmin(end, position))
  #gsub("[.]", "", s)
}, DATA[["germline_alignment"]], DATA[["v_germline_start"]], DATA[["v_germline_end"]])


# Replace N in sequence_v with the corresponding base from germline_v
DATA[["sequence_v2"]] <- mapply(function(seq, germline) {
  # Split into characters
  seq_chars <- unlist(strsplit(seq, ""))
  germ_chars <- unlist(strsplit(germline, ""))
  
  # Replace N with the corresponding germline base
  seq_chars[seq_chars == "N"] <- germ_chars[seq_chars == "N"]
  
  # Collapse back to a string 
  paste(seq_chars, collapse = "")
}, DATA[["sequence_v1"]], DATA[["germline_v"]])

DATA$v_call1 = DATA$v_call

DATA$v_call <- sapply(strsplit(DATA$v_call, ","), `[`, 1)

DATA$germline_v_extended <- mapply(function(seq_v, gene) {
  if (!gene %in% names(TRBV_GERM)) {
    # warning(paste("Missing gene in reference:", gene))
    return(seq_v)
  }
  
  ref <- as.character(TRBV_GERM[[gene]])
  return(ref)
}, DATA$germline_v, DATA$v_call)


DATA <- DATA[DATA[["v_mut"]] != 0,]
#DATA$v_call_gene <- sub("\\*.*", "", sub(",.*", "", DATA$v_call))

# Create summary table
#v_summary <- DATA %>%
#  group_by(sequence_v,germline_v, v_mut, v_call) %>%
#  summarise(count = n(), .groups = "drop")


# First prepare j_call_gene
DATA$j_call_gene <- sub("\\*.*", "", sub(",.*", "", DATA$j_call))

v_summary <- DATA %>%
  group_by(sequence_v2, v_call) %>%
  summarise(
    count = n(),
    unique_j_genes = n_distinct(j_call_gene),
    unique_cdr3_lengths = n_distinct(nchar(cdr3)),
    .groups = "drop"
  )%>%
  filter(!str_detect(v_call, ","))
  
# Then summarize
v_summary <- DATA %>%
  group_by(sequence_v2, germline_v_extended, v_call) %>%
  summarise(
    count = n(),
    unique_j_genes = n_distinct(j_call_gene),
    unique_cdr3_lengths = n_distinct(nchar(cdr3)),
    .groups = "drop"
  )%>%
  filter(!str_detect(v_call, ","))



#v_summary_top10 <- v_summary[order(-v_summary$count), ][1:10, ]

get_mutation_string <- function(seq_aln, germ_aln, start_pos) {

  if (is.na(seq_aln) || is.na(germ_aln)) {
    print("NA in one of them")
  }
  # Remove gaps for indexing
  seq_chars <- strsplit(seq_aln, "")[[1]]
  germ_chars <- strsplit(germ_aln, "")[[1]]
  
  mutation_list <- character()
  pos <- start_pos
  
  n <- min(length(seq_chars), length(germ_chars))
  
  for (i in seq_len(n)) {
    seq_char <- seq_chars[i]
    germ_char <- germ_chars[i]
    # If seq_char != germ_char, it's a mutation
    if (seq_char != germ_char) {
      mutation_list <- c(mutation_list, paste0(germ_char, pos, seq_char))
    }
    pos <- pos + 1
  }
  
  if (length(mutation_list) == 0) {
    return(NA)
  } else {
    return(paste(mutation_list, collapse = "_"))
  }
}

min_seq = 5 

v_summary_up <- v_summary[v_summary$count >= min_seq, ]

replace_dots_with_germline <- function(seq_aln, germ_aln) {
  if (is.na(seq_aln) || is.na(germ_aln)) {
    return(NA)
  }
  
  seq_chars <- strsplit(seq_aln, "")[[1]]
  germ_chars <- strsplit(germ_aln, "")[[1]]
  
  n <- min(length(seq_chars), length(germ_chars))
  # Replace "." with corresponding germline character
  seq_chars[seq_chars == "."] <- germ_chars[seq_chars == "."]
  
  return(paste(seq_chars, collapse = ""))
}

# Apply this before calling get_mutation_string
v_summary_up$sequence_v <- mapply(
  replace_dots_with_germline,
  v_summary_up$sequence_v2,
  v_summary_up$germline_v_extended
)


# Append matching reference bases to sequence_v2
v_summary_up$sequence_v_extended <- mapply(function(seq_v, gene) {
  if (!gene %in% names(TRBV_GERM)) {
    warning(paste("Missing gene in reference:", gene))
    return(seq_v)
  }
  ref <- as.character(TRBV_GERM[[gene]])
  ref_tail <- substr(ref, nchar(seq_v)+1, nchar(ref))
  extended_seq <- paste0(seq_v, ref_tail)

  return(extended_seq)
}, 
v_summary_up$sequence_v,
v_summary_up$v_call)


v_summary_up$mutation_string <- mapply(
  get_mutation_string,
  v_summary_up$sequence_v,
  v_summary_up$germline_v_extended,
  1
)

v_summary_up$novel_name <- ifelse(
  is.na(v_summary_up$mutation_string),
  v_summary_up$v_call,  # no mutations
  paste0(v_summary_up$v_call, "_", v_summary_up$mutation_string)
)

v_summary_up$novel_name <- as.character(v_summary_up$novel_name)

v_summary_up$novel_name

novel_igdiscover <- v_summary_up


novel_igdiscover <- novel_igdiscover[novel_igdiscover[["unique_cdr3_lengths"]] >= 2 & novel_igdiscover[["unique_j_genes"]] >= 2, ]

novel_igdiscover$NOTE <- NA
max_snp_position <- 316
rows2remove <- c() # initialize if not yet

for (i in seq_len(nrow(novel_igdiscover))) {
  allele <- novel_igdiscover$v_call[i]
  match_count <- length(grep(allele, data$v_call, fixed = TRUE))
  threshold <- match_count * 0.05
  if (novel_igdiscover$count[i] < threshold) {
    rows2remove <- c(rows2remove, i)
    novel_igdiscover$NOTE[i] <- paste0("Less than 5% (", round(threshold, 2), ") of the gene sequences were exact copies;")
  } else {
    novel_igdiscover$NOTE[i] <- paste0("More than 5% (", round(threshold, 2), ") of the gene sequences were exact copies;")
  }
  gene <- unlist(str_split(novel_igdiscover$novel_name[[i]], "[*]"))[1]
  gene <- unlist(str_split(gene, "-"))[1] # greps the family?
  ALLELES <- TRBV_GERM[grepl(gene, names(TRBV_GERM))]
  novel_imgt_seq <- novel_igdiscover$sequence_v[[i]]
  for (j in 1:length(ALLELES)) {
    if (grepl(ALLELES[[j]], novel_imgt_seq)) {
      new_name <- paste0(names(ALLELES[j]), "_", 
                         gsub(TRBV_GERM[names(ALLELES[j])], 
                              as.character(str_length(TRBV_GERM[names(ALLELES[j])])+1), 
                              novel_imgt_seq), 
                         str_length(novel_imgt_seq))
      name_index <- names(TRBV_GERM) == names(ALLELES[j])
      if (!grepl("_", names(ALLELES[j]))) {
        TRBV_GERM[names(ALLELES[j])] <- novel_imgt_seq
        names(TRBV_GERM)[name_index] <- new_name
      }
      rows2remove <- c(unlist(rows2remove), i)
      novel_igdiscover$NOTE[i] <- paste0(novel_igdiscover$NOTE[i],"Extantion of known allele;")
      next()
    }
    else {
      cutted_allele_seq <- substr(ALLELES[[j]], 1, max_snp_position)
      cutted_allele_seq <- gsub(".", "", cutted_allele_seq, fixed = T)
      cutted_allele_seq <- substr(cutted_allele_seq, 5, str_length(cutted_allele_seq))
      
      cutted_novel <- substr(novel_imgt_seq, 1, max_snp_position)
      cutted_novel <- gsub(".", "", cutted_novel, fixed = T)
      cutted_novel <- substr(cutted_novel, 5, str_length(cutted_novel))
      if(cutted_allele_seq == cutted_novel) {
        rows2remove <- c(unlist(rows2remove), i)
        novel_igdiscover$NOTE[i] <- paste0(novel_igdiscover$NOTE[i],"Identical to known allele;")
      }
    }
  }
}

if(length(rows2remove) > 0){
  novel_igdiscover <- novel_igdiscover[-rows2remove,]  
  }

write.table(novel_igdiscover, file = "igdiscover_novel_selected_igdiscover.tsv", sep = '\t', row.names = FALSE)

novel_genotype <- novel_igdiscover[, c("sequence_v_extended", "novel_name")]
#novel_genotype$novel_name <- gsub('-', '.', novel_genotype$novel_name, fixed = TRUE)
novel_genotype_vec <- setNames(novel_genotype$sequence_v_extended, novel_genotype$novel_name)

# Add only novel alleles not already in TRBV_GERM
for (novel_allele in names(novel_genotype_vec)) {
  if (!(novel_allele %in% names(TRBV_GERM))) {
    TRBV_GERM[novel_allele] <- novel_genotype_vec[novel_allele]
  }
}

personal_novel_fasta = paste0("","with_novel_V_ref.fasta")
tigger::writeFasta(TRBV_GERM, file = personal_novel_fasta)
