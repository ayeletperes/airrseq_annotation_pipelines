#!/usr/bin/env Rscript

log <- function(msg) {
  cat(paste(Sys.time(), ">", msg), sep = "\n")
}

log("Starting Undocumented Allele Inference")

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
v_germline_file <- args[2]
out_novel_file <- args[3]
out_novel_germline <- args[4]
cmd <- args[5]  # "library(tigger)" or source functions path
params <- args[6:length(args)]

log(paste("Loading environment:", cmd))
eval(parse(text = cmd))

suppressMessages(require(dplyr))

log("Reading input data and germline reference")
data <- data.table::fread(input_file, stringsAsFactors = FALSE, data.table = FALSE,
select = c('sequence_id', 'v_call', 'j_call', 'sequence_alignment', 'junction', 'junction_length', 'v_germline_end'))
vgerm <- tigger::readIgFasta(v_germline_file)

log("Transferring parameters")
num_threads   <- as.numeric(params[1])
germline_min  <- as.numeric(params[2])
min_seqs      <- as.numeric(params[3])
y_intercept   <- as.numeric(params[4])
alpha         <- as.numeric(params[5])
j_max         <- as.numeric(params[6])
min_frac      <- as.numeric(params[7])
auto_mutrange <- as.logical(params[8])
mut_range     <- as.numeric(unlist(strsplit(params[9], ":")))
mut_range     <- mut_range[1]:mut_range[2]
pos_range     <- as.numeric(unlist(strsplit(params[10], ":")))
pos_range     <- pos_range[1]:pos_range[2]

log("Running findNovelAlleles")
novel <- try(findNovelAlleles(
  data = data,
  germline_db = vgerm,
  v_call = 'v_call',
  j_call = 'j_call',
  seq = 'sequence_alignment',
  junction = 'junction',
  junction_length = 'junction_length',
  germline_min = germline_min,
  min_seqs = min_seqs,
  y_intercept = y_intercept,
  alpha = alpha,
  j_max = j_max,
  min_frac = min_frac,
  auto_mutrange = auto_mutrange,
  mut_range = mut_range,
  pos_range = pos_range,
  nproc = num_threads
))

if (!inherits(novel, 'try-error') && nrow(novel) > 0) {
  log("Filtering and selecting novel alleles")
  novel <- tigger::selectNovel(novel) %>%
    dplyr::distinct(novel_imgt, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(novel_imgt), nt_substitutions != '') %>%
    dplyr::mutate(gene = alakazam::getGene(germline_call, strip_d = FALSE)) %>%
    dplyr::group_by(gene) %>% dplyr::top_n(n = 2, wt = novel_imgt_count)

  log("Filtering padded alleles")
  Repeated_Read <- function(x, seq) {
    NT <- as.numeric(gsub('([0-9]+).*', '\\1', x))
    SNP <- gsub('.*>', '', x)
    OR_SNP <- gsub('[0-9]+([[:alpha:]]*).*', '\\1', x)
    seq <- c(substr(seq, (NT), (NT + 3)),
             substr(seq, (NT - 1), (NT + 2)),
             substr(seq, (NT - 2), (NT + 1)),
             substr(seq, (NT - 3), (NT)))
    PAT <- paste0(c(
      paste0(c(rep(SNP, 3), OR_SNP), collapse = ""),
      paste0(c(rep(SNP, 2), OR_SNP, SNP), collapse = ""),
      paste0(c(SNP, OR_SNP, rep(SNP, 2)), collapse = ""),
      paste0(c(OR_SNP, rep(SNP, 3)), collapse = "")
    ), collapse = '|')
    if (any(grepl(PAT, seq))) return(gsub(SNP, 'X', gsub(OR_SNP, 'z', seq[grepl(PAT, seq)])))
    else return(NA)
  }

  SNP_XXXX <- unlist(sapply(1:nrow(novel), function(i) {
    subs <- strsplit(novel[['nt_substitutions']][i], ',')[[1]]
    RR <- unlist(sapply(subs, Repeated_Read, seq = novel[['germline_imgt']][i], simplify = FALSE))
    RR <- RR[!is.na(RR)]
    length(RR) != 0
  }))
  novel <- novel[!SNP_XXXX, ]
  novel <- novel[!duplicated(novel[['polymorphism_call']]), ]

  if(nrow(novel) == 0) {
    log("No novel alleles found after filtering")
    quit(status = 0)
  }
  
  log(paste("Writing novel allele table:", out_novel_file))
  write.table(novel, file = out_novel_file, row.names = FALSE, quote = FALSE, sep = '\t')

  log(paste("Writing updated germline fasta: ", out_novel_germline, ".fasta\n", sep = ""))
  novel_v_germline <- setNames(gsub('-', '.', novel[['novel_imgt']], fixed = TRUE), novel[['polymorphism_call']])
  tigger::writeFasta(c(vgerm, novel_v_germline), paste0(out_novel_germline, '.fasta'))
} else {
  log("No novel alleles found or TIgGER failed")
  print(novel)
}
