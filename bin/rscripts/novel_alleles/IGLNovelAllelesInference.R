#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
v_germline_file <- args[2]
out_novel_file <- args[3]
out_novel_germline <- args[4]
cmd <- args[5]
params <- args[6:length(args)]

# Load Tigger or source environment
eval(parse(text = cmd))

suppressMessages(require(dplyr))

# Repeated read logic for light chain novel detection
Repeated_Read <- function(x, seq) {
  NT <- as.numeric(gsub('([0-9]+).*', '\\1', x))
  SNP <- gsub('.*>', '', x)
  OR_SNP <- gsub('[0-9]+([[:alpha:]]*).*', '\\1', x)
  seq <- c(substr(seq, (NT), (NT + 3)), substr(seq, (NT - 1), (NT + 2)), substr(seq, (NT - 2), (NT + 1)), substr(seq, (NT - 3), (NT)))
  PAT <- paste0(c(
    paste0(c(rep(SNP, 3), OR_SNP), collapse = ""),
    paste0(c(rep(SNP, 2), OR_SNP, SNP), collapse = ""),
    paste0(c(SNP, OR_SNP, rep(SNP, 2)), collapse = ""),
    paste0(c(OR_SNP, rep(SNP, 3)), collapse = "")
  ), collapse = '|')
  if (any(grepl(PAT, seq))) return(gsub(SNP, 'X', gsub(OR_SNP, 'z', seq[grepl(PAT, seq)])))
  else return(NA)
}

# Tigger upper bound helper
tigger_uper_bound <- function(DATA) {
  DATA <- DATA %>% dplyr::filter(!grepl(',', DATA[['v_call']])) %>%
    dplyr::select(v_call, v_germline_end) %>%
    dplyr::mutate(GENE = alakazam::getGene(v_call, strip_d = F)) %>%
    dplyr::group_by(GENE) %>% dplyr::mutate(RANGE = floor(quantile(v_germline_end, 0.95))) %>% dplyr::select(GENE, RANGE) %>%
    dplyr::slice(1) %>% dplyr::ungroup() %>% dplyr::group_by(RANGE) %>% dplyr::mutate(GENES = paste0(GENE, "[*]", collapse = "|")) %>%
    dplyr::slice(1) %>% dplyr::select(RANGE, GENES) %>% dplyr::ungroup()
  gene_range <- setNames(DATA[['GENES']], DATA[['RANGE']])
  return(gene_range)
}

# Load data and reference
data <- data.table::fread(input_file, stringsAsFactors = FALSE, data.table = FALSE,
select = c('sequence_id', 'v_call', 'j_call', 'sequence_alignment', 'junction', 'junction_length', 'v_germline_end'))
data <- data[substr(data[['sequence_alignment']], 1, 1) != ".",]
vgerm <- tigger::readIgFasta(v_germline_file)

# Parameters
num_threads <- as.numeric(params[1])
germline_min <- as.numeric(params[2])
min_seqs <- as.numeric(params[3])
y_intercept <- as.numeric(params[4])
alpha <- as.numeric(params[5])
j_max <- as.numeric(params[6])
min_frac <- as.numeric(params[7])
auto_mutrange <- as.logical(params[8])
mut_range <- as.numeric(unlist(strsplit(params[9], ":")))
mut_range <- mut_range[1]:mut_range[2]

# Get bounds
gene_range <- tigger_uper_bound(data)
novel <- data.frame()
for (i in seq_along(gene_range)) {
  upper_range <- as.numeric(names(gene_range)[i])
  genes <- gene_range[i]
  sub_ <- data[stringi::stri_detect_regex(data[["v_call"]], genes), ]
  if (nrow(sub_) != 0) {
    low_range <- min(sub_$v_germline_start)
    novel_df_tmp = try(findNovelAlleles(
      data = sub_,
      germline_db = vgerm,
      pos_range = low_range:upper_range,
      v_call = "v_call",
      j_call = "j_call",
      seq = "sequence_alignment",
      junction = "junction",
      junction_length = "junction_length",
      nproc = num_threads,
      germline_min = germline_min,
      min_seqs = min_seqs,
      y_intercept = y_intercept,
      alpha = alpha,
      j_max = j_max,
      min_frac = min_frac,
      auto_mutrange = auto_mutrange,
      mut_range = mut_range
    ))
    if (class(novel_df_tmp) != "try-error") {
      novel <- bind_rows(novel, novel_df_tmp)
    }
  }
}

if (!inherits(novel, 'try-error') && nrow(novel) > 0) {
  novel <- tigger::selectNovel(novel) %>%
    dplyr::distinct(novel_imgt, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(novel_imgt), nt_substitutions != '') %>%
    dplyr::mutate(gene = alakazam::getGene(germline_call, strip_d = FALSE)) %>%
    dplyr::group_by(gene) %>% dplyr::top_n(n = 2, wt = novel_imgt_count)
  SNP_XXXX <- unlist(sapply(seq_len(nrow(novel)), function(i) {
    subs <- strsplit(novel[['nt_substitutions']][i], ',')[[1]]
    RR <- unlist(sapply(subs, Repeated_Read, seq = novel[['germline_imgt']][i], simplify = FALSE))
    RR <- RR[!is.na(RR)]
    length(RR) != 0
  }))
  novel <- novel[!SNP_XXXX, ]
  novel <- novel[!duplicated(novel[['polymorphism_call']]), ]
  write.table(novel, file = out_novel_file, row.names = FALSE, quote = FALSE, sep = '\t')
  novel_v_germline <- setNames(gsub('-', '.', novel[['novel_imgt']], fixed = TRUE), novel[['polymorphism_call']])
  tigger::writeFasta(c(vgerm, novel_v_germline), paste0(out_novel_germline, '.fasta'))
}
