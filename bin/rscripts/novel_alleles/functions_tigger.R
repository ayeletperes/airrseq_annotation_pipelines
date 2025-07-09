
# The TIgGER Trifecta -----------------------------------------------------

#' Find novel alleles from repertoire sequencing data
#'
#' \code{findNovelAlleles} analyzes mutation patterns in sequences thought to
#' align to each germline allele in order to determine which positions
#' might be polymorphic.
#' 
#' The TIgGER allele-finding algorithm, briefly, works as follows:
#' Mutations are determined through comparison to the provided germline.
#' Mutation frequency at each *position* is determined as a function of
#' *sequence-wide* mutation counts. Polymorphic positions exhibit a high
#' mutation frequency despite sequence-wide mutation count. False positive of
#' potential novel alleles resulting from clonally-related sequences are guarded
#' against by ensuring that sequences perfectly matching the potential novel
#' allele utilize a wide range of combinations of J gene and junction length.
#' 
#' @param    data             \code{data.frame} containing repertoire data. See details.
#' @param    germline_db      vector of named nucleotide germline sequences
#'                            matching the V calls in \code{data}. These should be 
#'                            the gapped reference germlines used to make the V calls.
#' @param    v_call           name of the column in \code{data} with V allele calls. 
#'                            Default is \code{v_call}.    
#' @param    j_call           name of the column in \code{data} with J allele calls. 
#'                            Default is \code{j_call}. 
#' @param    seq              name of the column in \code{data} with the 
#'                            aligned, IMGT-numbered, V(D)J nucleotide sequence.
#'                            Default is \code{sequence_alignment}.
#' @param    junction         Junction region nucleotide sequence, which includes
#'                            the CDR3 and the two flanking conserved codons. Default
#'                            is \code{junction}.
#' @param    junction_length  Number of junction nucleotides in the junction sequence.
#'                            Default is \code{junction_length}.                        
#' @param    germline_min     the minimum number of sequences that must have a
#'                            particular germline allele call for the allele to
#'                            be analyzed.
#' @param    auto_mutrange    if \code{TRUE}, the algorithm will attempt to
#'                            determine the appropriate mutation range
#'                            automatically using the mutation count of the most
#'                            common sequence assigned to each allele analyzed.
#' @param    mut_range        range of mutations that samples may carry and
#'                            be considered by the algorithm.
#' @param    pos_range        range of IMGT-numbered positions that should be
#'                            considered by the algorithm.
#' @param    alpha            alpha value used for determining whether the 
#'                            fit y-intercept is greater than the \code{y_intercept}
#'                            threshold.
#' @param    y_intercept      y-intercept threshold above which positions should be
#'                            considered potentially polymorphic.
#' @param    j_max            maximum fraction of sequences perfectly aligning
#'                            to a potential novel allele that are allowed to
#'                            utilize to a particular combination of junction
#'                            length and J gene. The closer to 1, the less strict 
#'                            the filter for junction length and J gene diversity
#'                            will be.
#' @param    min_seqs         minimum number of total sequences (within the
#'                            desired mutational range and nucleotide range)
#'                            required for the samples to be considered.
#' @param    min_frac         minimum fraction of sequences that must have
#'                            usable nucleotides in a given position for that
#'                            position to considered.
#' @param    nproc            number of processors to use.
#'
#' @return
#' A \code{data.frame} with a row for each known allele analyzed.
#' Besides metadata on the the parameters used in the search, each row will have
#' either a note as to where the polymorphism-finding algorithm exited or a
#' nucleotide sequence for the predicted novel allele, along with columns providing
#' additional evidence.
#' 
#' The output contains the following columns:
#' \itemize{
#'   \item \code{germline_call}: The input (uncorrected) V call.
#'   \item \code{note}: Comments regarding the inferrence.
#'   \item \code{polymorphism_call}: The novel allele call.
#'   \item \code{nt_substitutions}: Mutations identified in the novel allele, relative
#'         to the reference germline (\code{germline_call})
#'   \item \code{novel_imgt}: The novel allele sequence.
#'   \item \code{novel_imgt_count}:  The number of times the sequence \code{novel_imgt} 
#'         is found in the input data. Considers the subsequence of \code{novel_imgt} 
#'         in the \code{pos_range}.
#'   \item \code{novel_imgt_unique_j}: Number of distinct J calls associated to \code{novel_imgt} 
#'         in the input data. Considers the subsequence of \code{novel_imgt} in the \code{pos_range}.       
#'   \item \code{novel_imgt_unique_cdr3}: Number of distinct CDR3 sequences associated
#'         with \code{novel_imgt} in the input data. Considers the subsequence of \code{novel_imgt} 
#'         in the \code{pos_range}.                                              
#'   \item \code{perfect_match_count}: Final number of sequences retained to call the new 
#'         allele. These are unique sequences that have V segments that perfectly match 
#'         the predicted germline in the \code{pos_range}.
#'   \item \code{perfect_match_freq}: \code{perfect_match_count / germline_call_count}
#'   \item \code{germline_call_count}: The number of sequences with the \code{germline_call} 
#'         in the input data that were initially considered for the analysis.
#'   \item \code{germline_call_freq}: The fraction of sequences with the \code{germline_call} 
#'         in the input data initially considered for the analysis.              
#'   \item \code{germline_imgt}: Germline sequence for \code{germline_call}.
#'   \item \code{germline_imgt_count}: The number of times the \code{germline_imgt} 
#'         sequence is found in the input data.
#'   \item \code{mut_min}: Minimum mutation considered by the algorithm.
#'   \item \code{mut_max}: Maximum mutation considered by the algorithm.
#'   \item \code{mut_pass_count}: Number of sequences in the mutation range.
#'   \item \code{pos_min}: First position of the sequence considered by the algorithm (IMGT numbering).
#'   \item \code{pos_max}: Last position of the sequence considered by the algorithm (IMGT numbering).
#'   \item \code{y_intercept}: The y-intercept above which positions were considered 
#'         potentially polymorphic.
#'   \item \code{y_intercept_pass}: Number of positions that pass the \code{y_intercept} threshold.
#'   \item \code{snp_pass}: Number of sequences that pass the \code{y_intercept} threshold and are
#'         within the desired nucleotide range (\code{min_seqs}).
#'   \item \code{unmutated_count}: Number of unmutated sequences.
#'   \item \code{unmutated_freq}: Number of unmutated sequences over \code{germline_imgt_count}.
#'   \item \code{unmutated_snp_j_gene_length_count}: Number of distinct combinations
#'         of SNP, J gene, and junction length.     
#'   \item \code{snp_min_seqs_j_max_pass}: Number of SNPs that pass both the \code{min_seqs} 
#'         and \code{j_max} thresholds.
#'   \item \code{alpha}: Significance threshold to be used when constructing the 
#'         confidence interval for the y-intercept.
#'   \item \code{min_seqs}: Input \code{min_seqs}. The minimum number of total sequences 
#'         (within the desired mutational range and nucleotide range) required 
#'         for the samples to be considered.
#'   \item \code{j_max}: Input \code{j_max}. The maximum fraction of sequences perfectly 
#'         aligning to a potential novel allele that are allowed to utilize to a particular 
#'         combination of junction length and J gene.
#'   \item \code{min_frac}: Input \code{min_frac}. The minimum fraction of sequences that must
#'         have usable nucleotides in a given position for that position to be considered.
#' }
#' 
#' The following comments can appear in the \code{note} column:
#' 
#' \itemize{
#'   \item \emph{Novel allele found}: A novel allele was detected.
#'   \item \emph{Plurality sequence too rare}: No sequence is frequent enough to pass 
#'         the J test (\code{j_max}).
#'   \item \emph{A J-junction combination is too prevalent}: Not enough J diversity (\code{j_max}).
#'   \item \emph{No positions pass y-intercept test}: No positions above \code{y_intercept}.
#'   \item \emph{Insufficient sequences in desired mutational range}: 
#'         \code{mut_range} and \code{pos_range}.
#'   \item \emph{Not enough sequences}: Not enough sequences in the desired mutational 
#'         range and nucleotide range (\code{min_seqs}).
#'   \item \emph{No unmutated versions of novel allele found}: All observed variants of the 
#'         allele are mutated.
#' }
#' 
#' @seealso \link{selectNovel} to filter the results to show only novel alleles.
#' \link{plotNovel} to visualize the data supporting any
#' novel alleles hypothesized to be present in the data and
#' \link{inferGenotype} and \link{inferGenotypeBayesian} to determine if the novel alleles are frequent
#' enought to be included in the subject's genotype.
#' 
#' @examples
#' \donttest{
#' # Note: In this example, with SampleGermlineIGHV,
#' # which contains reference germlines retrieved on August 2014,
#' # TIgGER finds the allele IGHV1-8*02_G234T. This allele
#' # was added to IMGT as IGHV1-8*03 on March 28, 2018.
#' 
#' # Find novel alleles and return relevant data
#' novel <- findNovelAlleles(AIRRDb, SampleGermlineIGHV)
#' selectNovel(novel)
#' }
#' 
#' @export
findNovelAlleles <- function(data, germline_db,
                             v_call="v_call",
                             j_call="j_call",
                             seq="sequence_alignment",
                             junction="junction",
                             junction_length="junction_length",
                             v_seq_length = "v_germline_end",
                             germline_min=200,
                             min_seqs=50,
                             auto_mutrange=TRUE,
                             mut_range=1:10,
                             pos_range=1:312,
                             y_intercept=0.125,
                             alpha=0.05,
                             j_max=0.15,
                             min_frac=0.75,
                             nproc=1) {
    . = idx = NULL
    suppressMessages(library(dplyr))
    suppressMessages(library(alakazam))
    suppressMessages(library(parallel))
    suppressMessages(library(foreach))    
    suppressMessages(library(stringi))

    # Keep only the db columns needed
    data <- data %>% 
        dplyr::select(!!!rlang::syms(c(seq,
                                       v_call,
                                       j_call,
                                       junction_length,
                                       junction,v_seq_length)))
    gc(verbose = FALSE)
    
    # Keep only the columns we need and clean up the sequences
    missing <- c(seq, v_call, j_call, junction_length,v_seq_length) %>%
        setdiff(colnames(data))
    if (length(missing) != 0) {
        stop("Could not find required columns in the input data:\n  ",
             paste(missing, collapse="\n  "))
    }
    empty_junctions <- sum(data[[junction_length]] == 0, na.rm=TRUE)
    if (empty_junctions > 0) {
        stop(empty_junctions, " sequences have junction ", "length of zero. ",
             "Please remove these sequences.")
    }
    germlines <- tigger::cleanSeqs(germline_db)
    names(germlines) <- getAllele(names(germlines), first=FALSE, strip_d=FALSE)
    data[[seq]] <- tigger::cleanSeqs(data[[seq]])
    
    
    # Find which rows' calls contain which germline alleles
    cutoff <-
        ifelse(germline_min < 1, round(nrow(data)*germline_min), germline_min)
    allele_groups <- sapply(names(germlines), grep, data[[v_call]], fixed=TRUE,
                           simplify=FALSE)
    names(allele_groups) <- names(germlines)
    allele_groups <- allele_groups[sapply(allele_groups, length) >= cutoff]
    if(length(allele_groups) == 0){
        stop_message <- paste("Not enough sample sequences were assigned to any germline:\n",
                              " (1) germline_min is too large or\n",
                              " (2) sequences names don't match germlines.")
        stop(stop_message)
    }
    allele_groups <- allele_groups[sortAlleles(names(allele_groups))]
    
    # Prepare for parallel processing
    nproc <- max(1,min(nproc, alakazam::cpuCount()-1))
    
    if(nproc == 1) {
        foreach::registerDoSEQ()
    } else {
        #cluster_type = ifelse(Sys.info()['sysname'] == "Windows", "PSOCK", "FORK")
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        parallel::clusterExport(cluster, list("allele_groups",
                                              "germlines",
                                              "v_call","j_call",
                                              "junction", "junction_length",
                                              "seq",
                                              "data",
                                              "min_seqs",
                                              "auto_mutrange",
                                              "mut_range",
                                              "pos_range",
                                              "y_intercept",
                                              "alpha",
                                              "j_max",
                                              "germline_min",
                                              "min_frac",
                                              "findLowerY",
                                              "mutationRangeSubset",
                                              "positionMutations",
                                              "superSubstring"),
                                envir=environment())
        doParallel::registerDoParallel(cluster)
    }
    
    out_list <- foreach(idx=iterators::icount(length(allele_groups))) %dopar% {
        library(dplyr)
        library(tigger)
        library(alakazam)
        # out_list <- lapply(1:length(allele_groups), function(idx) {  
        gc(verbose = FALSE) 
        #message(paste0("idx=",idx))
        # Subset of data being analyzed
        allele_name <- names(allele_groups)[idx]
        germline <- germlines[allele_name]
        indicies <- allele_groups[[allele_name]]
        db_subset <- data[indicies, ]
        
        # If mutrange is auto, find most popular mutation count and start from there
        gpm <- db_subset %>%
            dplyr::mutate(!!v_call := allele_name ) %>%
            tigger::getPopularMutationCount(germline,
                                    gene_min=0, seq_min=min_seqs,
                                    seq_p_of_max=1/8, full_return=TRUE,
                                    v_call=v_call,
                                    seq=seq)
        
        # Determine the mutation range(s) to scan
        mut_mins <- min(mut_range)
        if ( auto_mutrange & sum(gpm$mutation_count > 0) > 0 ){
            mut_mins <- c(mut_mins, gpm$mutation_count[gpm$mutation_count > 0]) %>%
                unique() %>%
                sort()
        }
        
        # Create the run's return object
        df_run_empty <- data.frame(germline_call = names(germline),
                                  note = "",
                                  polymorphism_call = NA,
                                  nt_substitutions=NA,
                                  novel_imgt = NA,
                                  novel_imgt_count=NA,
                                  novel_imgt_unique_j=NA,
                                  novel_imgt_unique_cdr3=NA,
                                  perfect_match_count = NA,
                                  perfect_match_freq = NA,                              
                                  germline_call_count = length(indicies),
                                  germline_call_freq = round(length(indicies)/nrow(data), 3),
                                  mut_min = NA,
                                  mut_max = NA,
                                  mut_pass_count=NA,
                                  germline_imgt = as.character(germline),
                                  germline_imgt_count=NA,
                                  pos_min = min(pos_range),
                                  pos_max = max(pos_range),
                                  y_intercept = y_intercept,
                                  y_intercept_pass = NA,
                                  snp_pass=NA,
                                  unmutated_count=NA,
                                  unmutated_freq=NA,
                                  unmutated_snp_j_gene_length_count=NA,
                                  snp_min_seqs_j_max_pass=NA,
                                  alpha = alpha,
                                  min_seqs = min_seqs,
                                  j_max = j_max,
                                  min_frac = min_frac,
                                  stringsAsFactors = FALSE)
        for (mut_min in rev(mut_mins)) {
            gc(verbose = FALSE)
            # message(paste0("|-- mut_min=",mut_min))
            if (mut_min == rev(mut_mins)[1]){
                df_run <- df_run_empty
            } else {
                df_run <- dplyr::bind_rows(df_run_empty, df_run)
            }
            mut_max <- mut_min + diff(range(mut_range))
            df_run$mut_min[1] <- mut_min
            df_run$mut_max[1] <- mut_max
            
            # If no sequence is frequent enough to pass the J test, give up now
            if(nrow(gpm) < 1) {
                df_run$note[1] <- "Plurality sequence too rare."
                if(mut_mins[1] == mut_min){
                    return(df_run)
                } else {
                    next
                }
            }
            
            # Add a mutation count column and filter out sequences not in our range
            db_subset_mm <- mutationRangeSubset(db_subset, germline,
                                               mut_min:mut_max, pos_range,
                                               seq=seq, v_seq_length = v_seq_length)
            df_run$mut_pass_count[1] <- nrow(db_subset_mm)
            
            if(nrow(db_subset_mm) < min_seqs){
                df_run$note[1] <- paste0("Insufficient sequences (",nrow(db_subset_mm),") in desired mutational range.")
                if(mut_mins[1] == mut_min){
                    return(df_run)
                } else {
                    next
                }
            }
            
            # Duplicate each sequence for all the positions to be analyzed
            # and find which positions are mutated
            pos_db <- positionMutations(db_subset_mm, germline, pos_range,
                                        seq=seq, v_seq_length = v_seq_length)
            
            # Find positional mut freq vs seq mut count
            pos_muts <- pos_db %>%
                dplyr::group_by(!!rlang::sym("POSITION")) %>%
                dplyr::mutate(PASS = mean(!! rlang::sym("OBSERVED")) >= min_frac) %>%
                dplyr::group_by(!!!rlang::syms(c("MUT_COUNT", "POSITION"))) %>%
                dplyr::summarise(POS_MUT_RATE = mean(!!rlang::sym("MUTATED"))*unique(!!rlang::sym("PASS")), count = n() ) %>% 
                dplyr::ungroup()   
            
            rm(pos_db)
            gc(verbose = FALSE)
            
            # Calculate y intercepts, find which pass the test
            pass_y <- pos_muts %>%
                dplyr::group_by(!!rlang::sym("POSITION")) %>%
                dplyr::summarise(Y_INT_MIN = findLowerY(!!rlang::sym("POS_MUT_RATE"),
                                                        !!rlang::sym("MUT_COUNT"),
                                                          mut_min, alpha, !!rlang::sym("count"))) %>%
                dplyr::filter(!!rlang::sym("Y_INT_MIN") > y_intercept)
            
            df_run$y_intercept_pass[1] <- nrow(pass_y)
            
            if(nrow(pass_y) < 1){
                df_run$note[1] <- "No positions pass y-intercept test."
                if(mut_mins[1] == mut_min){
                    return(df_run)
                } else {
                    next
                }
            }
            
            gl_substring <- superSubstring(germline, pass_y$POSITION)
            gl_minus_substring <- insertPolymorphisms(germline, pass_y$POSITION,
                                                     rep("N", nrow(pass_y)))
            
            # Find the potential SNP positions and remove anything that matches
            # the germline at all those positions or any combo that is too rare
            #message(idx)
            #message(pass_y$POSITION)
            db_y_subset_mm <- db_subset_mm %>%
                dplyr::group_by(1:n()) %>% dplyr::rowwise() %>%
                dplyr::mutate(SNP_STRING = superSubstring(!! rlang::sym(seq),
                                                            pass_y$POSITION), POSITION = !any(!!rlang::sym(v_seq_length)<pass_y$POSITION)) %>%
                dplyr::filter(!!rlang::sym("SNP_STRING") != gl_substring, !!rlang::sym("POSITION")) %>%
                dplyr::group_by(!!rlang::sym("SNP_STRING")) %>%
                dplyr::mutate(STRING_COUNT = n()) %>%
                dplyr::filter(!!rlang::sym("STRING_COUNT") >= min_seqs)
            
            df_run$snp_pass[1] <- nrow(db_y_subset_mm)
            
            if (nrow(db_y_subset_mm) < 1 ){
                df_run$note[1] <- paste("Position(s) passed y-intercept (",
                                       paste(pass_y$POSITION, collapse = ","),
                                       ") but the plurality sequence is too rare.",
                                       sep="")
                if(mut_mins[1] == mut_min){
                    return(df_run)
                } else {
                    next
                }
            }
            
            # Get mutation count at all positions that are not potential SNPs
            pads <- paste(rep("-", min(pos_range)-1), collapse="")
            db_y_subset_mm$MUT_COUNT_MINUS_SUBSTRING <-
                db_y_subset_mm %>% rowwise() %>%
                mutate(subseq = substring(!!rlang::sym(seq), min(pos_range), min(max(pos_range),!!rlang::sym(v_seq_length)))) %>%
                pull(subseq) %>%
                paste(pads, ., sep="") %>% 
                getMutatedPositions(gl_minus_substring) %>%
                sapply(length)
            
            # Keep only unmutated seqences and then find the counts of J and
            # junction length for each of the SNP strings, and then check to
            # see which pass the j/junction and count requirements
            db_y_summary0 <- db_y_subset_mm %>%
                dplyr::filter(!!rlang::sym("MUT_COUNT_MINUS_SUBSTRING") == 0)
            
            df_run$unmutated_count[1] <- nrow(db_y_summary0)
            
            db_y_summary0 <- db_y_summary0 %>%
                dplyr::mutate(J_GENE = getGene(!! rlang::sym(j_call))) %>%
                dplyr::group_by(!!! rlang::syms(c("SNP_STRING", "J_GENE", junction_length))) %>%
                dplyr::summarise(COUNT = n())
            
            df_run$unmutated_snp_j_gene_length_count[1] <- nrow(db_y_summary0)
            
            db_y_summary0 <- db_y_summary0 %>%
                dplyr::group_by(!!rlang::sym("SNP_STRING")) %>%
                dplyr::mutate(FRACTION = !! rlang::sym("COUNT")/sum(!! rlang::sym("COUNT"))) %>%
                dplyr::summarise(TOTAL_COUNT = sum(!! rlang::sym("COUNT")), MAX_FRAC = max(!! rlang::sym("FRACTION")))
            
            if(nrow(db_y_summary0) < 1){
                df_run$note[1] <- paste("Position(s) passed y-intercept (",
                                       paste(pass_y$POSITION, collapse = ","),
                                       ") but no unmutated versions of novel allele",
                                       " found.", sep="")
                if(mut_mins[1] == mut_min){
                    return(df_run)
                } else {
                    next
                }
            }
            
            # db_y_summary = db_y_summary0 %>%
            #   filter_(~TOTAL_COUNT >= min_seqs & MAX_FRAC <= j_max)
            
            min_seqs_pass <- db_y_summary0$TOTAL_COUNT >= min_seqs
            j_max_pass <- db_y_summary0$MAX_FRAC <= j_max
            
            db_y_summary <- db_y_summary0[min_seqs_pass & j_max_pass, , drop=FALSE]
            
            df_run$snp_min_seqs_j_max_pass[1] <- nrow(db_y_summary)
            
            if(nrow(db_y_summary) < 1){
                msg <- c(NA, NA)
                names(msg) <- c("j_max", "min_seqs")
                
                if (sum(min_seqs_pass) == 0) {
                    msg['min_seqs'] <- paste0("Not enough sequences (maximum total count is ",
                                              max(db_y_summary0$TOTAL_COUNT),
                                              ").")
                }
                
                if (sum(j_max_pass) == 0) {
                    msg['j_max'] <- paste0("A J-junction combination is too prevalent (",
                                           round(100*max(db_y_summary0$MAX_FRAC),1),"% of sequences).")
                }
                
                msg <- paste(na.omit(msg), collapse=" and ")
                df_run$note[1] <- paste("Position(s) passed y-intercept (",
                                       paste(pass_y$POSITION, collapse = ","),
                                       ") but ",
                                       msg,".", sep="")
                df_run$perfect_match_count[1] <- max(db_y_summary0$TOTAL_COUNT)
                df_run$perfect_match_freq[1] <- df_run$perfect_match_count[1]/df_run$germline_call_count[1]
                if(mut_mins[1] == mut_min){
                    return(df_run)
                } else {
                    next
                }
            }
            
            germ_nts <- unlist(strsplit(gl_substring,""))
            for (r in 1:nrow(db_y_summary)) {
                if (r > 1){
                    df_run <- dplyr::bind_rows(df_run[1,], df_run)
                }
                # Create the new germline
                snp_nts <- unlist(strsplit(db_y_summary$SNP_STRING[r],""))
                remain_mut <- db_y_summary$SNP_STRING[r] %>%
                    getMutatedPositions(gl_substring) %>%
                    unlist() %>%
                    unique()
                germ <- insertPolymorphisms(germline, pass_y$POSITION, snp_nts)
                is_known_allele <- germ == germlines
                if (sum(is_known_allele) == 0 ) {
                    names(germ) <- mapply(paste, germ_nts[remain_mut],
                                         pass_y$POSITION[remain_mut],
                                         snp_nts[remain_mut], sep="") %>%
                        paste(collapse="_") %>%
                        paste(names(germline), ., sep="_")
                } else {
                    # If the match is with duplicated sequences in the reference germlines,
                    # use the first
                    known_allele_names <- sortAlleles(names(germlines)[is_known_allele],
                                                      method="position")
                    names(germ) = known_allele_names[1]
                }
                # Save the new germline to our data frame               
                df_run$polymorphism_call[1] <- names(germ)
                df_run$novel_imgt[1] <-  as.character(germ)
                df_run$perfect_match_count[1] <- db_y_summary$TOTAL_COUNT[r]
                df_run$perfect_match_freq[1] <- df_run$perfect_match_count[1]/df_run$germline_call_count[1]
                df_run$note[1] = "Novel allele found!"
            }
            
        } # end for each starting mutation counts
        return(df_run)
        
    } # end foreach allele
    
        if(nproc > 1) { stopCluster(cluster) }
    out_df <- dplyr::bind_rows(out_list)
    getMuSpec <- function(poly_call, germ_call) {
        sapply(1:length(poly_call), function(i){
            p <- gsub(germ_call[i], "", poly_call[i], fixed = T)
            p <- strsplit(p,"_")[[1]][-1]
            m <- gsub("([[:alpha:]])([[:digit:]]*)([[:alpha:]])", "\\2\\1>\\3", p)
            paste(m, collapse=",")
        })
    }
    
    # The number of records in the sequence dataset matching 
    # each exact novel_imgt sequence
    getDbMatch <- function(novel_imgt) {
        novel_imgt <- stri_sub(novel_imgt,  min(pos_range), max(pos_range))
        novel_imgt <- stri_replace_all_regex(novel_imgt, "[-\\.]","")
        data_seq <- stri_sub(data[[seq]],  min(pos_range), max(pos_range))
        data_seq <- stri_replace_all_regex(data_seq,"[-\\.]","")
        sapply(novel_imgt, function(n) {
            # sum(grepl(n, data_seq))
            sum(stri_detect_fixed(data_seq,n))
        })
    }
    
    # The number of distinct J in the sequence dataset associated 
    # with the exact novel_imgt sequence
    getNumJ <- function(novel_imgt) {
        novel_imgt <- stri_sub(novel_imgt,  min(pos_range), max(pos_range))
        novel_imgt <- stri_replace_all_regex(novel_imgt, "[-\\.]","")
        data_seq <- stri_sub(data[[seq]],  min(pos_range), max(pos_range))
        data_seq <- stri_replace_all_regex(data_seq,"[-\\.]","")      
        sapply(novel_imgt, function(n) {
            imgt_idx <- stri_detect_fixed(data_seq,n)
            length(unique(getGene(data[[j_call]][imgt_idx])))
        })
    }
    
    
    # The number of distinct CDR3 in the sequence dataset associated 
    # with the exact novel_imgt sequence
    getNumCDR3 <- function(novel_imgt) {
        novel_imgt <- stri_sub(novel_imgt,  min(pos_range), max(pos_range))
        novel_imgt <- stri_replace_all_regex(novel_imgt, "[-\\.]","")
        data_seq <- stri_sub(data[[seq]],  min(pos_range), max(pos_range))
        data_seq <- stri_replace_all_regex(data_seq,"[-\\.]","")            
        sapply(novel_imgt, function(n) {
            imgt_idx <- stri_detect_fixed(data_seq,n)
            seq <- data[[junction]][imgt_idx]
            seq <- stri_sub(seq, 4, stringi::stri_length(seq) - 3)
            length(unique(seq))
        })
    }
    
    idx <- which(!is.na(out_df$novel_imgt))
    if (length(idx)>0) {
        out_df$nt_substitutions[idx] <- getMuSpec(out_df$polymorphism_call[idx],
                                                  out_df$germline_call[idx])
        out_df$novel_imgt_count[idx] <- getDbMatch(out_df$novel_imgt[idx])
        out_df$novel_imgt_unique_j[idx] <- getNumJ(out_df$novel_imgt[idx])
        if (junction %in% colnames(data)) {
            out_df$novel_imgt_unique_cdr3[idx] <- getNumCDR3(out_df$novel_imgt[idx])
        }
    }
    out_df$germline_imgt_count <- getDbMatch(out_df$germline_imgt)
    out_df$unmutated_freq <- out_df$unmutated_count/out_df$germline_call_count

    return(out_df)
}

#' Select rows containing novel alleles
#' 
#' \code{selectNovel} takes the result from \link{findNovelAlleles} and
#' selects only the rows containing unique, novel alleles.
#' 
#' @details
#' If, for instance, subject has in his genome \code{IGHV1-2*02} and a novel 
#' allele equally close to \code{IGHV1-2*02} and \code{IGHV1-2*05}, the novel allele may be
#' detected by analyzing sequences that best align to either of these alleles.
#' If \code{keep_alleles} is \code{TRUE}, both polymorphic allele calls will
#' be retained. In the case that multiple mutation ranges are checked for the
#' same allele, only one mutation range will be kept in the output.
#' 
#' @param   novel           \code{data.frame} of the type returned by
#'                          \link{findNovelAlleles}.
#' @param   keep_alleles    \code{logical} indicating if different alleles
#'                          leading to the same novel sequence should be kept.
#'                          See Details.
#'                        
#' @return  A \code{data.frame} containing only unique, novel alleles (if any)
#' that were in the input.
#' 
#' @examples
#' novel <- selectNovel(SampleNovel)
#' 
#' @export
selectNovel <- function(novel, keep_alleles=FALSE) {
    # Remove non-novel rows
    novel <- filter(novel, !is.na(!!rlang::sym("novel_imgt")))
    
    if (keep_alleles) {
        novel < novel %>% 
            group_by(!!rlang::sym("germline_call"))
    }
    novel_set <- novel %>%
        distinct(!!rlang::sym("novel_imgt"), .keep_all=TRUE) %>%
        ungroup()
    
    return(novel_set)
}

#' Visualize evidence of novel V alleles
#'
#' \code{plotNovel} is be used to visualize the evidence of any novel V
#' alleles found using \link{findNovelAlleles}. It can also be used to
#' visualize the results for alleles that did
#' 
#' @details
#' The first panel in the plot shows, for all sequences which align to a particular 
#' germline allele, the mutation frequency at each postion along the aligned 
#' sequence as a function of the sequence-wide mutation count. Each line is a position.
#' Positions that contain polymorphisms (rather than somatic hypermutations) 
#' will exhibit a high apparent mutation frequency for a range of 
#' sequence-wide mutation counts. The positions are color coded as follows:
#' 
#' \itemize{
#'   \item  red:    the position(s) pass(ess) the novel allele test 
#'   \item  yellow: the position(s) pass(ess) the y-intercept test but not
#'                         other tests
#'   \item  blue:   the position(s) didn't pass the y-intercept test and 
#'                         was(were) not further considered
#' }
#'  
#' The second panel shows the nucleotide usage at each of the polymorphic positions
#' as a function of sequence-wide mutation count. If no polymorphisms were identified,
#' the panel will show the mutation count.
#' 
#' To avoid cases where a clonal expansion might lead to a false positive, TIgGER examines
#' the combinations of J gene and junction length among sequences which perfectly 
#' match the proposed germline allele. Clonally related sequences usually share 
#' the same V gene, J gene and junction length. Requiring the novel allele
#' to be found in different combinations of J gene and junction lengths
#' is a proxy for requiring it to be found in different clonal lineages.
#' 
#' @param    data             \code{data.frame} containing repertoire data. See
#'                            \link{findNovelAlleles} for details.
#' @param    novel_row        single row from a data frame as output by
#'                            \link{findNovelAlleles} that contains a
#'                            polymorphism-containing germline allele.
#' @param    v_call           name of the column in \code{data} with V allele
#'                            calls. Default is \code{v_call}.
#' @param    j_call           name of the column in \code{data} with J allele calls. 
#'                            Default is \code{j_call}. 
#' @param    seq              name of the column in \code{data} with the 
#'                            aligned, IMGT-numbered, V(D)J nucleotide sequence.
#'                            Default is \code{sequence_alignment}.
#' @param    junction         Junction region nucleotide sequence, which includes
#'                            the CDR3 and the two flanking conserved codons. Default
#'                            is \code{junction}.
#' @param    junction_length  number of junction nucleotides in the junction sequence.
#'                            Default is \code{junction_length}.                        
#' @param    ncol             number of columns to use when laying out the plots.
#' @param    multiplot        whether to return one single plot (\code{TRUE}) or a list 
#'                            with the three individual plots (\code{FALSE}).
#' @examples
#' # Plot the evidence for the first (and only) novel allele in the example data
#' novel <- selectNovel(SampleNovel)
#' plotNovel(AIRRDb, novel[1, ], v_call="v_call", j_call="j_call", 
#'           seq="sequence_alignment", junction="junction", junction_length="junction_length", 
#'           multiplot=TRUE)
#' 
#' @export
plotNovel <- function(data, novel_row, v_call="v_call", j_call="j_call",
                      seq="sequence_alignment",
                      junction="junction", junction_length="junction_length",v_seq_length="v_germline_end",
                      ncol=1, multiplot=TRUE) {
    . = NULL
    
    # Use the data frame
    if(length(novel_row) > 0) {
        if(is.data.frame(novel_row) & nrow(novel_row) == 1) {
            pos_range <- novel_row$pos_min:novel_row$pos_max
            germline <- novel_row$germline_imgt
            names(germline) <- novel_row$germline_call
            mut_range <- novel_row$mut_min[1]:novel_row$mut_max[1]
            novel_imgt <- novel_row$novel_imgt
            names(novel_imgt) <- novel_row$polymorphism_call
            min_frac <- novel_row$min_frac
            note <- novel_row$note
        } else {
            stop("novel_row is not a data frame with only one row.")
        }
    }
    
    germline <- cleanSeqs(germline)
    data[[seq]] <- cleanSeqs(data[[seq]])
    
    # Extract sequences assigned to the germline, determine which
    # have an appropriate range of mutations, and find the mutation
    # frequency of each position
    db_subset <- data %>%
        select(!!!rlang::syms(c(seq, v_call, j_call, junction_length))) %>%
        filter(grepl(names(germline),  .data[[v_call]], fixed=TRUE))
    pos_db <- db_subset %>%  
        mutationRangeSubset(germline, mut_range, pos_range, seq=seq, v_seq_length = v_seq_length)
    if (nrow(pos_db) == 0) {
        warning(paste0("Insufficient sequences (",nrow(pos_db),") in desired mutational range."))
        return (invisible(NULL))
    }
    pos_db <- pos_db %>%
        positionMutations(germline, pos_range,seq=seq, v_seq_length=v_seq_length)
    pos_muts <- pos_db %>%
        group_by(!!rlang::sym("POSITION")) %>%
        mutate(PASS = mean(!!rlang::sym("OBSERVED")) >= min_frac) %>%
        group_by(!!!rlang::syms(c("MUT_COUNT", "POSITION"))) %>%
        summarise(POS_MUT_RATE = mean(!!rlang::sym("MUTATED"))*unique(!!rlang::sym("PASS"))) %>% 
        ungroup()
    
    # Label the polymorphic positions as such
    pass_y <- unlist(strsplit(names(novel_imgt), "_"))[-1] %>%
        gsub("[^0-9]", "", .) %>%
        as.numeric()
    p_y_f <- unlist(strsplit(names(novel_imgt), "_"))[-1] %>%
        gsub("[0-9]+.", "", .)
    p_y_t <- unlist(strsplit(names(novel_imgt), "_"))[-1] %>%
        gsub(".[0-9]+", "", .)
    # Parse the note to find positions that passed y intercept if no novel found
    if (length(pass_y) == 0 & grepl("Position\\(s\\) passed y-intercept", note)) {
        pass_y <- note %>% 
            gsub("Position\\(s\\) passed y-intercept \\(", "", .) %>%
            gsub("\\).*", "", .) %>% strsplit(",") %>% unlist %>% as.numeric
        p_y_f <- sapply(pass_y, function (x) substring(germline, x, x))
        p_y_t <- gsub(".", "?", p_y_f)
    }
    
    to_from <- paste(paste("Position", pass_y), paste(paste(p_y_f, "->"), p_y_t))
    names(to_from) <- pass_y
    pos_muts <- pos_muts %>%
        mutate(Polymorphic = ifelse(!!rlang::sym("POSITION") %in% pass_y, "True", "False"))
    
    pads <- paste(rep("-", min(pos_range)-1), collapse="")
    db_subset$MUT_COUNT_NOVEL <- db_subset[[seq]] %>%
        substring(min(pos_range), max(pos_range)) %>%
        paste(pads, ., sep="") %>%
        getMutatedPositions(novel_imgt) %>%
        sapply(length)
    db_subset <- db_subset %>%
        filter(!!rlang::sym("MUT_COUNT_NOVEL") == 0) %>%
        mutate(J_GENE = getGene(!!rlang::sym(j_call)))
    if (nrow(db_subset) == 0) {
        warning(paste0("Insufficient sequences (",nrow(db_subset),") with MUT_COUNT_NOVEL == 0."))
        return (invisible(NULL))
    }
    db_subset[[junction_length]] <- db_subset[[junction_length]] %>%
        factor(levels=min(db_subset[[junction_length]]):max(db_subset[[junction_length]]))
    pos_muts$Polymorphic <- pos_muts$Polymorphic %>%
        factor(levels = c("False", "True"))
    pos_db$NT <- pos_db$NT %>%
        factor(levels = names(DNA_COLORS))
    pos_muts$GERMLINE <- names(germline)
    
    # MAKE THE FIRST PLOT
    if (!is.na(novel_imgt)) {
        POLYCOLORS <- setNames(DNA_COLORS[c(4,3)], c("False", "True")) # blue #3C88EE, red #EB413C
        p1 <- ggplot(pos_muts, aes_string(x="MUT_COUNT", y="POS_MUT_RATE", 
                                          group="POSITION", color="Polymorphic")) +
            geom_line(data=filter(pos_muts, !!rlang::sym("Polymorphic") == "False"), size=0.75) +
            geom_line(data=filter(pos_muts, !!rlang::sym("Polymorphic") == "True"), size=0.75) +
            facet_grid(GERMLINE ~ .) +
            scale_color_manual(values = POLYCOLORS) +
            ylim(0,1) +
            xlab("Mutation Count (Sequence)") +
            ylab("Mutation Frequency (Position)") +
            theme_bw() +
            theme(legend.position=c(0.5,0.9), legend.justification=c(0.5,1),
                  legend.background=element_rect(fill = "transparent")) +
            guides(color = guide_legend(ncol = 2, reverse = TRUE))
    } else{
        POLYCOLORS <- setNames(DNA_COLORS[c(4,2)], c("False", "True")) # blue #3C88EE, yellow #FFB340
        p1 <- ggplot(pos_muts, aes_string(x="MUT_COUNT", y="POS_MUT_RATE", 
                                          group="POSITION", color="Polymorphic")) +
            geom_line(size=0.75) +
            facet_grid(GERMLINE ~ .) +
            scale_color_manual(values=POLYCOLORS) +
            ylim(0, 1) +
            xlab("Mutation Count (Sequence)") +
            ylab("Mutation Frequency (Position)") +
            theme_bw() +
            theme(legend.position=c(0.5, 0.9), legend.justification=c(0.5, 1),
                  legend.background=element_rect(fill="transparent")) +
            guides(color = guide_legend("Passed y-intercept test",
                                        ncol = 2, reverse = TRUE))
    }
    # MAKE THE SECOND PLOT
    p2_data <- mutate(filter(pos_db, !!rlang::sym("POSITION") %in% pass_y),
                      POSITION = to_from[as.character(!!rlang::sym("POSITION"))])
    positions <- unique(p2_data$POSITION)
    numeric_positions <- as.numeric(sub("Position ([0-9]+) .+","\\1",positions))
    
    p2_data$POSITION <- factor(p2_data$POSITION,
                               levels=positions[order(numeric_positions)],
                               ordered = TRUE)
    
    if (nrow(p2_data)) {
        p2 <- ggplot(p2_data, aes_string(x="MUT_COUNT", fill="NT")) +
            geom_bar(width=0.9) +
            guides(fill = guide_legend("Nucleotide", ncol=4)) +
            xlab("Mutation Count (Sequence)") + 
            ylab("Sequence Count") +
            scale_fill_manual(values=DNA_COLORS, breaks=names(DNA_COLORS),
                              drop=FALSE) +
            theme_bw() +
            theme(legend.position=c(1,1), legend.justification=c(1,1),
                  legend.background=element_rect(fill="transparent"))
    } else {
        p2_data <- mutate(filter(pos_db,
                                  !!rlang::sym("POSITION") %in% names(which.max(table(pos_db$POSITION)))),
                          POSITION = "No positions pass y-intercept test.")
        p2 <- ggplot(p2_data, aes_string(x="MUT_COUNT")) +
            geom_bar(width=0.9) +
            xlab("Mutation Count (Sequence)") + ylab("Sequence Count") +
            theme_bw() +
            theme(legend.position=c(1,1), legend.justification=c(1,1),
                  legend.background=element_rect(fill = "transparent"))
    }
    
    if (length(positions) < 5) {
        p2 <- p2 +
            facet_grid(POSITION ~ .) 
        ncols <- 1
        nrows <- length(positions)
    } else {
        ncols <- min(length(positions),6)
        nrows <- ceiling(length(positions)/ncols)
        p2 <- p2 +
            facet_wrap(POSITION ~ ., ncol=ncols)
    }
        
    # MAKE THE THIRD PLOT
    p3 <- ggplot(db_subset, aes_string(x=junction_length, fill="J_GENE")) +
        geom_bar(width=0.9) +
        guides(fill=guide_legend("J Gene", ncol=2)) +
        xlab("Junction Length") + 
        ylab("Unmutated Sequence Count") +
        theme_bw() +
        theme(legend.position=c(1, 1), legend.justification=c(1, 1),
              legend.background=element_rect(fill="transparent"))
    
    p2_height <- max(1,0.6*nrows )

    heights <- c(1, p2_height, 1)
    if (multiplot) {
        multiplot(p1, p2, p3, cols = ncol, heights=heights)         
    } else {
        list(p1, p2, p3)
    }
}

#' Infer a subject-specific genotype using a frequency method
#'
#' \code{inferGenotype} infers an subject's genotype using a frequency method.
#' The genotype is inferred by finding the minimum number set of alleles that 
#' can explain the majority of each gene's calls. The most common allele of 
#' each gene is included in the genotype first, and the next most common allele 
#' is added until the desired fraction of alleles can be explained. In this 
#' way, mistaken allele calls (resulting from sequences which
#' by chance have been mutated to look like another allele) can be removed.
#' 
#' @details
#' Allele calls representing cases where multiple alleles have been
#' assigned to a single sample sequence are rare among unmutated
#' sequences but may result if nucleotides for certain positions are
#' not available. Calls containing multiple alleles are treated as
#' belonging to all groups. If \code{novel} is provided, all
#' sequences that are assigned to the same starting allele as any
#' novel germline allele will have the novel germline allele appended
#' to their assignent prior to searching for unmutated sequences.
#'           
#' @param    data                 \code{data.frame} containing V allele
#'                                calls from a single subject.
#' @param    germline_db          named vector of sequences containing the
#'                                germline sequences named in
#'                                \code{allele_calls}. Only required if
#'                                \code{find_unmutated} is \code{TRUE}.
#' @param    novel                optional \code{data.frame} of the type
#'                                novel returned by
#'                                \link{findNovelAlleles} containing
#'                                germline sequences that will be utilized if
#'                                \code{find_unmutated} is \code{TRUE}. See
#'                                Details.
#' @param    v_call               column in \code{data} with V allele calls.
#'                                Default is \code{"v_call"}.                            
#' @param    seq                  name of the column in \code{data} with the 
#'                                aligned, IMGT-numbered, V(D)J nucleotide sequence.
#'                                Default is \code{sequence_alignment}.
#' @param    fraction_to_explain  the portion of each gene that must be
#'                                explained by the alleles that will be included
#'                                in the genotype.
#' @param    gene_cutoff          either a number of sequences or a fraction of
#'                                the length of \code{allele_calls} denoting the
#'                                minimum number of times a gene must be
#'                                observed in \code{allele_calls} to be included
#'                                in the genotype.
#' @param    find_unmutated       if \code{TRUE}, use \code{germline_db} to
#'                                find which samples are unmutated. Not needed
#'                                if \code{allele_calls} only represent
#'                                unmutated samples.
#' 
#' @return
#' A \code{data.frame} of alleles denoting the genotype of the subject containing 
#' the following columns:
#'           
#' \itemize{
#'   \item \code{gene}: The gene name without allele.
#'   \item \code{alleles}: Comma separated list of alleles for the given \code{gene}.
#'   \item \code{counts}: Comma separated list of observed sequences for each 
#'         corresponding allele in the \code{alleles} list.
#'   \item \code{total}: The total count of observed sequences for the given \code{gene}.
#'   \item \code{note}: Any comments on the inferrence.
#' }
#'           
#' @note
#' This method works best with data derived from blood, where a large
#' portion of sequences are expected to be unmutated. Ideally, there
#' should be hundreds of allele calls per gene in the input.
#' 
#' @seealso \link{plotGenotype} for a colorful visualization and
#'          \link{genotypeFasta} to convert the genotype to nucleotide sequences.
#'          See \link{inferGenotypeBayesian} to infer a subject-specific genotype 
#'          using a Bayesian approach.
#' 
#' @examples
#' # Infer IGHV genotype, using only unmutated sequences, including novel alleles
#' inferGenotype(AIRRDb, germline_db=SampleGermlineIGHV, novel=SampleNovel,
#'               find_unmutated=TRUE)
#' 
#' @export
inferGenotype <- function(data, germline_db=NA, novel=NA, v_call="v_call", 
                          seq="sequence_alignment",
                          fraction_to_explain=0.875, gene_cutoff=1e-4, 
                          find_unmutated=TRUE) {
    
    . = NULL
    allele_calls = getAllele(data[[v_call]], first=FALSE, strip_d=FALSE)
    # Find the unmutated subset, if requested
    if (find_unmutated) {
        if(is.na(germline_db[1])){
            stop("germline_db needed if find_unmutated is TRUE")
        }
        if (!is.null(nrow(novel))) {
            novel <- filter(novel, !is.na(!!rlang::sym("polymorphism_call"))) %>%
                select(!!!rlang::syms(c("germline_call", "polymorphism_call", "novel_imgt")))
            if (nrow(novel) > 0) {
                # Extract novel alleles if any and add them to germline_db
                novel_gl <- novel$novel_imgt
                names(novel_gl) <- novel$polymorphism_call
                germline_db <- c(germline_db, novel_gl)
                # Add the novel allele calls to allele calls of the same starting allele
                for(r in 1:nrow(novel)){
                    ind <- grep(novel$germline_call[r], allele_calls, fixed=TRUE)
                    allele_calls[ind] <- allele_calls[ind] %>%
                        sapply(paste, novel$polymorphism_call[r], sep=",")
                }
            }
        }
        # Find unmutated sequences
        allele_calls <- findUnmutatedCalls(allele_calls,
                                          as.character(data[[seq]]),
                                          germline_db)
        if(length(allele_calls) == 0){
            stop("No unmutated sequences found! Set 'find_unmutated' to 'FALSE'.")
        }
    }
    
    # Find which rows' calls contain which genes
    cutoff <- ifelse(gene_cutoff < 1, length(allele_calls)*gene_cutoff, gene_cutoff)
    gene_regex <- allele_calls %>% strsplit(",") %>% unlist() %>%
        getGene(strip_d=FALSE) %>%  unique() %>% paste("\\*", sep="")
    gene_groups <- sapply(gene_regex, grep, allele_calls, simplify=FALSE)
    names(gene_groups) <- gsub("\\*", "", gene_regex, fixed=TRUE)
    gene_groups <- gene_groups[sapply(gene_groups, length) >= cutoff]
    gene_groups <- gene_groups[sortAlleles(names(gene_groups))]
    
    # Make a table to store the resulting genotype
    gene <- names(gene_groups)
    alleles <- counts <- note <- rep("", length(gene))
    total <- sapply(gene_groups, length)
    genotype <- cbind(gene, alleles, counts, total, note)
    
    # For each gene, find which alleles to include
    for (g in gene) {
        # Keep only the part of the allele calls that uses the gene being analyzed
        ac <- allele_calls[gene_groups[[g]]] %>%
            strsplit(",") %>%
            lapply(function(x) x[grep(paste(g, "\\*", sep=""), x)]) %>%
            sapply(paste, collapse=",")
        target <- ceiling(fraction_to_explain*length(ac)) # how many we need to explain
        t_ac <- table(ac) # table of allele calls
        potentials <- unique(unlist(strsplit(names(t_ac),","))) # potential alleles
        # One allele? Easy!
        if (length(potentials) == 1 | length(t_ac) == 1) {
            genotype[genotype[,"gene"]==g,"alleles"] <- gsub("[^d\\*]*[d\\*]","",potentials )[1]
            genotype[genotype[,"gene"]==g,"counts"] <- t_ac
        } else {
            # More alleles? Let's find the fewest that can explain the needed fraction
            # Make a table of which alleles can explain which calls
            regexpotentials <- paste(gsub("\\*","\\\\*", potentials),"$",sep="")
            regexpotentials <- 
                paste(regexpotentials,gsub("\\$",",",regexpotentials),sep="|")
            tmat = 
                sapply(regexpotentials, function(x) grepl(x, names(t_ac),fixed=FALSE))
            seqs_expl = as.data.frame(apply(tmat, 2, function(x) x*t_ac))
            colnames(seqs_expl) = potentials
            
            # Cycle through the table, including alleles to explain more sequences,
            # until we explain enough sequences
            included <- counts <- character(0)
            tot_expl = 0
            while(tot_expl < target){
                allele_tot <- apply(seqs_expl, 2, sum)
                included <- c(included, names(which.max(allele_tot)))
                counts <- c(counts, max(allele_tot))
                tot_expl <- max(allele_tot)  + tot_expl
                seqs_expl <- seqs_expl[which(seqs_expl[,which.max(allele_tot)]==0),]
            }
            genotype[genotype[,"gene"]==g,"alleles"] <-
                paste(gsub("[^d\\*]*[d\\*]","",included ),collapse=",")
            genotype[genotype[,"gene"]==g,"counts"] <-
                paste(counts,collapse=",")
        }
    }
    
    geno <- as.data.frame(genotype, stringsAsFactors = FALSE)
    
    # Check for indistinguishable calls
    if (find_unmutated == TRUE) {
        seqs <- genotypeFasta(geno, germline_db)
        dist_mat <- seqs %>%
            sapply(function(x) sapply((getMutatedPositions(seqs, x)), length)) %>%
            as.matrix
        rownames(dist_mat) <- colnames(dist_mat)
        for (i in 1:nrow(dist_mat)){ dist_mat[i,i] = NA }
        same <- which(dist_mat == 0, arr.ind=TRUE)
        if (nrow(same) > 0 ) {
            for (r in 1:nrow(same)) {
                inds <- as.vector(same[r,])
                geno[getGene(rownames(dist_mat)[inds][1]),]$note <-
                    paste(rownames(dist_mat)[inds], collapse=" and ") %>%
                    paste("Cannot distinguish", .)
            }
        }
    }
    rownames(geno) <- NULL
    
    return(geno)
}


#' Show a colorful representation of a genotype
#'
#' \code{plotGenotype} plots a genotype table.
#' 
#' @param    genotype     \code{data.frame} of alleles denoting a genotype, 
#'                        as returned by \link{inferGenotype}.
#' @param    facet_by     column name in \code{genotype} to facet the plot by. 
#'                        if \code{NULL}, then do not facet the plot. 
#' @param    gene_sort    string defining the method to use when sorting alleles.
#'                        if \code{"name"} then sort in lexicographic order. If
#'                        \code{"position"} then sort by position in the locus, as
#'                        determined by the final two numbers in the gene name.
#' @param    text_size    point size of the plotted text.
#' @param    silent       if \code{TRUE} do not draw the plot and just return the ggplot
#'                        object; if \code{FALSE} draw the plot.
#' @param    ...          additional arguments to pass to ggplot2::theme.
#' 
#' @return  A ggplot object defining the plot.
#' 
#' @seealso \link{inferGenotype}
#' 
#' @examples
#' # Plot genotype
#' plotGenotype(SampleGenotype)
#' 
#' # Facet by subject
#' genotype_a <- genotype_b <- SampleGenotype
#' genotype_a$SUBJECT <- "A"
#' genotype_b$SUBJECT <- "B"
#' geno_sub <- rbind(genotype_a, genotype_b)
#' plotGenotype(geno_sub, facet_by="SUBJECT", gene_sort="pos")
#' 
#' @export
plotGenotype <- function(genotype, facet_by=NULL, gene_sort=c("name", "position"), 
                         text_size=12, silent=FALSE, ...) {
    # Check arguments
    gene_sort <- match.arg(gene_sort)
    
    # Split genes' alleles into their own rows
    alleles = strsplit(genotype$alleles, ",")
    geno2 = genotype
    r = 1
    for (g in 1:nrow(genotype)){
        for(a in 1:length(alleles[[g]])) {
            geno2[r, ] = genotype[g, ]
            geno2[r, ]$alleles = alleles[[g]][a]
            r = r + 1
        }
    }
    
    # Set the gene order
    geno2$gene = factor(geno2$gene, 
                        levels=rev(sortAlleles(unique(geno2$gene), method=gene_sort)))
    
    # Create the base plot
    p = ggplot(geno2, aes_string(x="gene", fill="alleles")) +
        theme_bw() +
        theme(axis.ticks=element_blank(),
              axis.text.x=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              text=element_text(size=text_size),
              strip.background=element_blank(),
              strip.text=element_text(face="bold")) +
        geom_bar(position="fill") +
        coord_flip() + xlab("Gene") + ylab("") +
        scale_fill_hue(name="Allele", h=c(0, 270), h.start=10)
    
    # Plot, with facets by SUBJECT if that column is present
    if (!is.null(facet_by)) {
        p = p + facet_grid(paste0(".~", facet_by))
    }
    
    # Add additional theme elements
    p = p + do.call(theme, list(...))
    
    # Plot
    if (!silent) { plot(p) }
    
    invisible(p)
}

#' Return the nucleotide sequences of a genotype
#'
#' \code{genotypeFasta} converts a genotype table into a vector of nucleotide
#' sequences.
#' 
#' @param    genotype     \code{data.frame} of alleles denoting a genotype, 
#'                        as returned by \link{inferGenotype}.
#' @param    germline_db  vector of named nucleotide germline sequences
#'                        matching the alleles detailed in \code{genotype}.
#' @param    novel        an optional \code{data.frame} containing putative
#'                        novel alleeles of the type returned by
#'                        \link{findNovelAlleles}.
#' 
#' @return   A named vector of strings containing the germline nucleotide
#'           sequences of the alleles in the provided genotype.
#' 
#' @seealso \link{inferGenotype}
#' 
#' @examples
#' # Find the sequences that correspond to the genotype
#' genotype_db <- genotypeFasta(SampleGenotype, SampleGermlineIGHV, SampleNovel)
#' 
#' @export
genotypeFasta <- function(genotype, germline_db, novel=NA){
    if(!is.null(nrow(novel))){
        # Extract novel alleles if any and add them to germline_db
        novel <- filter(novel, !is.na(!!rlang::sym("polymorphism_call"))) %>%
            select(!!!rlang::syms(c("germline_call", "polymorphism_call", "novel_imgt")))
        if(nrow(novel) > 0){
            novel_gl <- novel$novel_imgt
            names(novel_gl) <- novel$polymorphism_call
            germline_db <- c(germline_db, novel_gl)
        }
    }

    genotype$gene <- gsub("[Dd]\\*","*",genotype$gene)
    g_names <- names(germline_db)
    names(g_names) <- gsub("[Dd]\\*", "*", names(germline_db))
    table_calls <- mapply(paste, genotype$gene, strsplit(genotype$alleles, ","),
                         sep="*")
    seqs <- germline_db[as.vector(g_names[unlist(table_calls)])]
    if(sum(is.na(seqs)) > 0){
        stop("The following genotype alleles were not found in germline_db: ",
             paste(unlist(table_calls)[which(is.na(seqs))], collapse = ", "))
    }
    return(seqs)
}

#' Correct allele calls based on a personalized genotype
#'
#' \code{reassignAlleles} uses a subject-specific genotype to correct
#' correct preliminary allele assignments of a set of sequences derived
#' from a single subject.
#' 
#' @details
#' In order to save time, initial gene assignments are preserved and
#' the allele calls are chosen from among those provided in \code{genotype_db},
#' based on a simple alignment to the sample sequence.
#' 
#' @param    data          \code{data.frame} containing V allele calls from a
#'                         single subject and the sample IMGT-gapped V(D)J sequences under
#'                         \code{seq}.
#' @param    genotype_db   vector of named nucleotide germline sequences
#'                         matching the calls detailed in \code{allele_calls}
#'                         and personalized to the subject
#' @param    v_call        name of the column in \code{data} with V allele
#'                         calls. Default is \code{v_call}.    
#' @param    seq           name of the column in \code{data} with the 
#'                         aligned, IMGT-numbered, V(D)J nucleotide sequence.
#'                         Default is SEQUENCE_IMGT                                      
#' @param    method        method to use when realigning sequences to
#'                         the genotype_db sequences. Currently, only \code{"hammming"}
#'                         (for Hamming distance) is implemented.
#' @param    path          directory containing the tool used in the
#'                         realignment method, if needed. Hamming distance does
#'                         not require a path to a tool.
#' @param    keep_gene     string indicating if the gene (\code{"gene"}), 
#'                         family (\code{"family"}) or complete repertoire
#'                         (\code{"repertoire"}) assignments should be performed. 
#'                         Use of \code{"gene"} increases speed by minimizing required number of 
#'                         alignments, as gene level assignments will be maintained when possible.
#' 
#' @return   A modifed input \code{data.frame} containing the best allele call from 
#'           among the sequences listed in \code{genotype_db} in the 
#'           \code{v_call_genotyped} column.
#' 
#' @examples
#' # Extract the database sequences that correspond to the genotype
#' genotype_db <- genotypeFasta(SampleGenotype, SampleGermlineIGHV, novel=SampleNovel)
#' 
#' # Use the personlized genotype to determine corrected allele assignments
#' output_db <- reassignAlleles(AIRRDb, genotype_db, v_call="v_call",
#'                              seq="sequence_alignment")
#' 
#' @export
reassignAlleles <- function(data, genotype_db, v_call="v_call",
                            seq="sequence_alignment",
                            method="hamming", path=NA,
                            keep_gene=c("gene", "family", "repertoire")){
    # Check arguments    
    keep_gene <- match.arg(keep_gene)
    
    # Extract data subset and prepare output vector
    v_sequences = as.character(data[[seq]])
    v_calls = getAllele(data[[v_call]], first=FALSE, strip_d=FALSE)
    v_call_genotyped = rep("", length(v_calls))
    
    if (keep_gene == "gene") { 
        v = getGene(v_calls, first = TRUE, strip_d=FALSE)
        geno = getGene(names(genotype_db),strip_d=TRUE)
        names(geno) = names(genotype_db)
    } else if (keep_gene == "family") {
        v <- getFamily(v_calls, first = TRUE, strip_d = FALSE)
        geno = getFamily(names(genotype_db),strip_d=TRUE)
        names(geno) = names(genotype_db)
    } else if (keep_gene == "repertoire") {
        v <- rep(v_call, length(v_calls))
        geno = rep(v_call, length(genotype_db))
        names(geno) = names(genotype_db)      
    } else {
        stop("Unknown keep_gene value: ", keep_gene)
    }
    
    # keep_gene == FALSE
    # Find which genotype genes/families are homozygous and assign those alleles first
    hetero = unique(geno[which(duplicated(geno))])
    homo = geno[!(geno %in% hetero)]
    homo_alleles = names(homo)
    names(homo_alleles) = homo
    homo_calls_i = which(v %in% homo)
    v_call_genotyped[homo_calls_i] = homo_alleles[v[homo_calls_i]]
    
    # Now realign the heterozygote sequences to each allele of that gene
    for (het in hetero){
        ind = which(v %in% het)
        if (length(ind) > 0){
            het_alleles = names(geno[which(geno == het)])
            het_seqs = genotype_db[het_alleles]
            if(method == "hamming"){
                dists = lapply(het_seqs, function(x)
                    sapply(getMutatedPositions(v_sequences[ind], x, match_instead=FALSE),
                           length))
                dist_mat = matrix(unlist(dists), ncol = length(het_seqs))
            } else {
                stop("Only Hamming distance is currently supported as a method.")
            }
            # The sapply-apply approach could become problematic when nrow(dist_mat)
            # is 1 and min(best_match) has multiple values, due to the fact that R 
            # does not always keep data structures unmutable
            # Explicitly specifying a list and subsequently keeping it as a list by
            # using lapply avoids that problem
            best_match = vector("list", length=nrow(dist_mat))
            for (i in 1:nrow(dist_mat)) {
                best_match[[i]] = which(dist_mat[i, ]==min(dist_mat[i, ]))
            }
            best_alleles = lapply(best_match, function(x) het_alleles[x])
            v_call_genotyped[ind] = unlist(lapply(best_alleles, paste, collapse=","))
        }
    }
    
    # Now realign the gene-not-in-genotype calls to every genotype allele
    hetero_calls_i = which(v %in% hetero)
    not_called = setdiff(1:length(v), c(homo_calls_i, hetero_calls_i))
    if(length(not_called)>1){
        if(method ==  "hamming"){
            dists = lapply(genotype_db, function(x)
                sapply(getMutatedPositions(v_sequences[not_called], x, match_instead=FALSE),
                       length))
            dist_mat = matrix(unlist(dists), ncol = length(genotype_db))
        } else {
            stop("Only Hamming distance is currently supported as a method.")
        }
        # The sapply-apply approach could become problematic when nrow(dist_mat)
        # is 1 and min(best_match) has multiple values, due to the fact that R 
        # does not always keep data structures unmutable
        # Explicitly specifying a list and subsequently keeping it as a list by
        # using lapply avoids that problem
        best_match = vector("list", length=nrow(dist_mat))
        for (i in 1:nrow(dist_mat)) {
            best_match[[i]] = which(dist_mat[i, ]==min(dist_mat[i, ]))
        }
        best_alleles = lapply(best_match, function(x) names(genotype_db[x]))
        v_call_genotyped[not_called] = unlist(lapply(best_alleles, paste, collapse=","))
    }
    
    if (all(v_call_genotyped == data[[v_call]])) {
        msg <- ("No allele assignment corrections made.") 
        if (all(v %in% homo) & length(hetero) > 0) {
            keep_opt <- eval(formals(reassignAlleles)$keep_gene)
            i <- match(keep_gene, keep_opt)
            rec_opt <- paste(keep_opt[(i+1):length(keep_opt)], collapse = ", ")
            msg <- paste(msg, "Consider setting keep_gene to one of:", rec_opt)
        }
        warning(msg)
    }
    
    data$v_call_genotyped <- v_call_genotyped
    
    return(data)
}


# Other Mutation-Related Functions ----------------------------------------

#' Find the location of mutations in a sequence
#'
#' \code{getMutatedPositions} takes two vectors of aligned sequences and
#' compares pairs of sequences. It returns a list of the nucleotide positions of
#' any differences.
#' 
#' @param    samples        vector of strings respresenting aligned sequences
#' @param    germlines      vector of strings respresenting aligned sequences
#'                          to which \code{samples} will be compared. If only
#'                          one string is submitted, it will be used for all
#'                          \code{samples}.
#' @param    ignored_regex  regular expression indicating what characters
#'                          should be ignored (such as gaps and N nucleotides).
#' @param    match_instead  if \code{TRUE}, the function returns the positions
#'                          that are the same instead of those that are
#'                          different.
#' @return   A list of the nucleotide positions of any differences between the
#'           input vectors.
#' 
#' @examples
#' # Create strings to act as a sample sequences and a reference sequence
#' seqs <- c("----GATA", "GAGAGAGA", "TANA")
#' ref <- "GATAGATA"
#' 
#' # Find the differences between the two
#' getMutatedPositions(seqs, ref)
#' 
#' @export
getMutatedPositions <- function(samples, germlines, ignored_regex="[\\.N-]",
                                match_instead=FALSE) {
    
    # If only one germline sequence is given, use it for all the sample seqs
    if(length(germlines) == 1){ germlines = rep(germlines, length(samples)) }
    if(length(samples) != length(germlines)) {
        stop("Number of input sequences does not match number of germlines.")
    }
    
    # Truncate each pair of sequences to the length of the shorter
    germ_mins = lapply(germlines, nchar)
    samp_mins = lapply(samples, nchar)
    min_lens = mapply(min, germ_mins, samp_mins)
    germ = toupper(mapply(substr, germlines, 1, min_lens, SIMPLIFY=FALSE))
    samp = toupper(mapply(substr, samples, 1, min_lens, SIMPLIFY=FALSE))
    
    # Calculate poisitions of mutations (or matches), ignoring gaps, Ns, and CDR3
    samp_char = strsplit(samp,"")
    germ_char = strsplit(germ,"")
    if(!match_instead){
        muts = lapply(mapply("!=", samp_char, germ_char, SIMPLIFY=FALSE), which)
    } else {
        muts = lapply(mapply("==", samp_char, germ_char, SIMPLIFY=FALSE), which)
    }
    ignore_germ = gregexpr(ignored_regex, germ)
    ignore_samp = gregexpr(ignored_regex, samp)
    ignore = mapply(c, ignore_germ, ignore_samp, SIMPLIFY=FALSE)
    
    muts = mapply(function(x, y) x[!x%in%y], muts, ignore, SIMPLIFY=FALSE)
    return(muts)
}


#' Determine the mutation counts from allele calls
#'
#' \code{getMutCount} takes a set of nucleotide sequences and their allele calls
#' and determines the distance between that seqeunce and any germline alleles
#' contained within the call
#' 
#' @param    samples       vector of IMGT-gapped sample V sequences
#' @param    allele_calls  vector of strings respresenting Ig allele calls for
#'                         the sequences in \code{samples}, where multiple
#'                         calls are separated by a comma
#' @param    germline_db   vector of named nucleotide germline sequences
#'                         matching the calls detailed in \code{allele_calls}
#' 
#' @return   A list equal in length to \code{samples}, containing the Hamming
#'           distance to each germline allele contained within each call within
#'           each element of \code{samples}
#' 
#' @examples
#' # Insert a mutation into a germline sequence
#' s2 <- s3 <- SampleGermlineIGHV[1]
#' stringi::stri_sub(s2, 103, 103) <- "G"
#' stringi::stri_sub(s3, 107, 107) <- "C"
#' 
#' sample_seqs <- c(SampleGermlineIGHV[2], s2, s3)
#' 
#' # Pretend that one sample sequence has received an ambiguous allele call
#' sample_alleles <- c(paste(names(SampleGermlineIGHV[1:2]), collapse=","),
#'                     names(SampleGermlineIGHV[2]),
#'                     names(SampleGermlineIGHV[1]))
#' 
#' # Compare each sequence to its assigned germline(s) to determine the distance
#' getMutCount(sample_seqs, sample_alleles, SampleGermlineIGHV)
#' 
#' @export
getMutCount <- function(samples, allele_calls, germline_db){
    
    call_list = strsplit(allele_calls, ",")
    
    germline_list = lapply(call_list, function(x) germline_db[x])
    
    mut_pos_list = list()
    mut_count_list = list()
    # First, find mutations of all sequences with call count of 1
    call_count = sapply(germline_list, length)
    cc1 = which(call_count == 1)
    if (length(cc1) > 0) {
        mut_pos_list[cc1] = getMutatedPositions(samples[cc1],
                                                unlist(germline_list[cc1]))
        mut_count_list[cc1] = lapply(mut_pos_list[cc1], length)
    }
    # Then find mutations of all sequences with call count > 1
    ccm = which(call_count > 1)
    if (length(ccm) > 0){
        mut_pos_list[ccm] = mapply(getMutatedPositions,
                                   germline_list[ccm], samples[ccm],
                                   SIMPLIFY=FALSE)
        mut_count_list[ccm] = lapply(mut_pos_list[ccm],
                                     function(x) lapply(x,length))
    }
    
    return(mut_count_list)
}

#' Determine which calls represent an unmutated allele
#'
#' \code{findUnmutatedCalls} determines which allele calls would represent a 
#' perfect match with the germline sequence, given a vector of allele calls and
#' mutation counts. In the case of multiple alleles being assigned to a
#' sequence, only the subset that would represent a perfect match is returned.
#' 
#' @param    allele_calls   vector of strings respresenting Ig allele calls,
#'                          where multiple calls are separated by a comma.
#' @param    germline_db    vector of named nucleotide germline sequences
#' @param    sample_seqs    V(D)J-rearranged sample sequences matching the order
#'                          of the given \code{allele_calls}.
#' 
#' @return   A vector of strings containing the members of \code{allele_calls}
#'           that represent unmutated sequences.
#' 
#' @examples
#' # Find which of the sample alleles are unmutated
#' calls <- findUnmutatedCalls(AIRRDb$v_call, AIRRDb$sequence_alignment, 
#'                             germline_db=SampleGermlineIGHV)
#' 
#' @export
findUnmutatedCalls <- function(allele_calls, sample_seqs, germline_db){
    . = NULL
    allele_calls = getAllele(allele_calls, first = FALSE)
    sample_seqs = as.character(sample_seqs)
    
    # Remove calls not in germline_db
    not_in_db = allele_calls %>%
        strsplit(",") %>%
        unlist %>%
        setdiff(names(germline_db))
    no_call = which(allele_calls == "") 
    in_db = not_in_db %>%
        sapply(grep, allele_calls, fixed=TRUE) %>%
        unlist() %>%
        c(no_call) %>%
        unique() %>%
        setdiff(1:length(allele_calls), .)
    allele_calls = allele_calls[in_db]
    sample_seqs = sample_seqs[in_db]
    
    mut_counts = getMutCount(sample_seqs, allele_calls, germline_db)
    
    # Find which seqs are unmutated and which of the allele calls that represents
    unmut_i = which(sapply(mut_counts, function(x) min(unlist(x))) == 0)
    which_no_muts = sapply(mut_counts, function(x) grep("^0$", unlist(x)) )
    unmut_alleles = rep("", length(allele_calls))
    
    # How many alleles represent perfect matches?
    n_gl_unmut = sapply(which_no_muts, length)
    
    one_unmut = which(n_gl_unmut == 1)
    split_names = strsplit(allele_calls, ",")
    if (length(one_unmut) > 0){
        inds = unlist(which_no_muts[one_unmut])
        unmut_alleles[one_unmut] = mapply("[", split_names[one_unmut], inds)
    }
    
    more_unmut = which(n_gl_unmut > 1)
    if (length(more_unmut) > 0){
        inds = which_no_muts[more_unmut]      
        unmut_multi = mapply(function(x,y) x[unlist(y)], split_names[more_unmut],
                             inds, SIMPLIFY = FALSE)
        unmut_alleles[more_unmut] = sapply(unmut_multi, paste, collapse=",")  
    }
    
    unmut_alleles = unmut_alleles[unmut_i]
    
    return(unmut_alleles)
    
}

#' Find mutation counts for frequency sequences
#'
#' \code{getPopularMutationCount} determines which sequences occur frequently
#' for each V gene and returns the mutation count of those sequences.
#' 
#' @param  data          \code{data.frame} in the Change-O format. See
#'                       \link{findNovelAlleles} for a list of required
#'                       columns.
#' @param  germline_db   named list of IMGT-gapped germline sequences.
#' @param  v_call        name of the column in \code{data} with V allele calls. 
#'                       Default is \code{v_call}.    
#' @param  seq           name of the column in \code{data} with the 
#'                       aligned, IMGT-numbered, V(D)J nucleotide sequence.
#'                       Default is \code{sequence_alignment}.
#' @param  gene_min      portion of all unique sequences a gene must
#'                       constitute to avoid exclusion.
#' @param  seq_min       number of copies of the V that must be present for
#'                       to avoid exclusion.
#' @param  seq_p_of_max  ror each gene, the fraction of the most common V sequence
#'                       count that a sequence must meet to avoid exclusion.
#' @param  full_return   if \code{TRUE}, will return all \code{data} columns and
#'                       will include sequences with mutation count < 1.
#' 
#' @return  A data frame of genes that have a frequent sequence mutation count
#'          above 1.
#' 
#' @seealso \link{getMutatedPositions} can be used to find which positions
#'          of a set of sequences are mutated.
#' 
#' @examples
#' getPopularMutationCount(AIRRDb, SampleGermlineIGHV)
#' 
#' @export
getPopularMutationCount <- function(data, germline_db, 
                                    v_call="v_call",
                                    seq="sequence_alignment",
                                    gene_min = 1e-03,
                                    seq_min = 50, seq_p_of_max = 1/8,
                                    full_return = FALSE){
    modified_db <- data %>%
        mutate(v_gene = getGene(!!rlang::sym(v_call))) %>%
        group_by(!!rlang::sym("v_gene")) %>%
        mutate(v_gene_n = n()) %>%
        group_by(1:n()) %>%
        mutate(v_sequence_imgt = substring(!!rlang::sym(seq), 1, 312)) %>%
        # Count occurence of each unique IMGT-gapped V sequence
        group_by(!!!rlang::syms(c("v_gene", "v_sequence_imgt"))) %>%
        mutate(v_sequence_imgt_n = n()) %>%
        # Determine count of most common sequence
        group_by(!!rlang::sym("v_gene")) %>%
        mutate(v_sequence_imgt_n_max = max(!!rlang::sym("v_sequence_imgt_n"))) %>%
        # Remove rare V genes, rare sequences, and sequences not making up a
        # sufficient proportion of sequences as compared to the most common
        ungroup %>%
        distinct(!!rlang::sym("v_sequence_imgt"), .keep_all = TRUE) %>%
        filter(!!rlang::sym("v_gene_n") >= (nrow(data)*gene_min)) %>%
        filter(!!rlang::sym("v_sequence_imgt_n") >= seq_min) %>%
        mutate(v_sequence_imgt_p_max = !!rlang::sym("v_sequence_imgt_n")/!!rlang::sym("v_sequence_imgt_n_max")) %>%
        filter(!!rlang::sym("v_sequence_imgt_p_max") >= seq_p_of_max)
    # Determine the mutation counts of the V sequences and append them to the db
    mutation_count <- getMutCount(modified_db$v_sequence_imgt,
                                 modified_db[[v_call]],
                                 germline_db) %>% 
        sapply(function(x) min(unlist(x)))
    if (length(mutation_count)==0){
        mutation_count <- integer(0)
    }
    merged_db <- bind_cols(modified_db, data.frame(mutation_count))
    # Strip down the data frame before returning it
    if (!full_return) {
        merged_db <- merged_db %>%
            filter(mutation_count > 0) %>%
            select(!!!rlang::syms(c("v_gene", "mutation_count")))
    }
    return(merged_db)
}

#' Insert polymorphisms into a nucleotide sequence
#'
#' \code{insertPolymorphisms} replaces nucleotides in the desired locations of a
#' provided sequence.
#' 
#' @param    sequence     starting nucletide sequence.
#' @param    positions    numeric vector of positions which to be changed.
#' @param    nucleotides  character vector of nucletides to which to change the
#'                        positions.
#'                        
#' @return   A sequence with the desired nucleotides in the provided locations.
#' 
#' @examples
#' insertPolymorphisms("HUGGED", c(1, 6, 2), c("T", "R", "I")) 
#' 
#' @export
insertPolymorphisms <- function(sequence, positions, nucleotides) {
    
    if(length(positions) != length(nucleotides)){
        stop("Number of nucleotides and number of positions do not match.")
    }
    names(positions) = nucleotides
    for (i in 1:length(positions)){
        substr(sequence, positions[i], positions[i]) = names(positions[i])
    }
    
    return(sequence)
}

# Formatting and Cleanup --------------------------------------------------

#' Read immunoglobulin sequences
#'
#' \code{readIgFasta} reads a fasta-formatted file of immunoglobulin (Ig)
#' sequences and returns a named vector of those sequences.
#' 
#' @param    fasta_file       fasta-formatted file of immunoglobuling sequences.
#' @param    strip_down_name  if \code{TRUE}, will extract only the allele name
#'                            from the strings fasta file's sequence names.
#' @param    force_caps       if \code{TRUE}, will force nucleotides to
#'                            uppercase.
#'                            
#' @return   Named vector of strings respresenting Ig alleles.
#' 
#' @seealso  \link{writeFasta} to do the inverse.
#' 
#' @examples 
#' \dontrun{
#' # germlines <- readIgFasta("ighv.fasta")
#' }
#' 
#' @export
readIgFasta <- function(fasta_file, strip_down_name=TRUE, force_caps=TRUE) {
    all_char = readChar(fasta_file, file.info(fasta_file)$size)
    split_by_sequence = strsplit(all_char, "[ \t\r\n\v\f]?>")
    add_name_break = sapply(split_by_sequence, function(x) sub("[\r\n]",">",x))
    cleaned_up = sapply(add_name_break, function(x) gsub("[ \t\r\n\v\f]", "", x))
    broken_names = sapply(cleaned_up, strsplit, ">")
    
    seqs = sapply(broken_names, "[", 2)
    seq_names = sapply(broken_names, "[", 1)
    if(force_caps) { seqs = toupper(seqs) }
    if(strip_down_name){ seq_names = getAllele(seq_names, strip_d=FALSE) }
    names(seqs) = seq_names
    
    return(seqs[which(!is.na(seqs))])
}

#' Write to a fasta file
#'
#' \code{writeFasta} writes a named vector of sequences to a file in fasta
#' format.
#' 
#' @param    named_sequences  vector of named string representing sequences
#' @param    file             the name of the output file.
#' @param    width            the number of characters to be printed per line.
#'                            if not between 1 and 255, width with be infinite.
#' @param    append           \code{logical} indicating if the output should be
#'                            appended to \code{file} instead of overwriting it
#' 
#' @return   A named vector of strings respresenting Ig alleles.
#' 
#' @seealso  \link{readIgFasta} to do the inverse.
#' 
#' @examples
#' \dontrun{
#' # writeFasta(germlines, "ighv.fasta")
#' }
#' 
#' @export
writeFasta <- function(named_sequences, file, width=60, append=FALSE){
    . = NULL
    seq_names = names(named_sequences) %>%
        paste(">", ., "\n", sep="")
    seqs = as.character(named_sequences)
    if(is.numeric(width) & width > 0 & width < 256){
        width_regex = paste("(.{", width, ",", width, "})", sep="")
        seqs = gsub(width_regex, "\\1\n", seqs)
    }
    seqs = seqs %>%
        paste("\n", sep="") %>%
        gsub("\n\n", "\n", .)
    paste(seq_names, seqs, sep="", collapse="") %>%
        cat(file=file, append=append)
} 

#' Update IGHV allele names
#'
#' \code{updateAlleleNames} takes a set of IGHV allele calls and replaces any
#' outdated names (e.g. IGHV1-f) with the new IMGT names.
#' 
#' @param    allele_calls  vector of strings respresenting IGHV allele names.
#' 
#' @return   Vector of strings respresenting updated IGHV allele names.
#' 
#' @note
#' IGMT has removed \code{IGHV2-5*10} and \code{IGHV2-5*07} as it has determined they
#' are actually alleles \code{02} and \code{04}, respectively. The updated allele 
#' names are based on IMGT release 2014-08-4.
#' 
#' @references
#' \enumerate{
#'   \item Xochelli et al. (2014) Immunoglobulin heavy variable (IGHV) genes
#'         and alleles: new entities, new names and implications for research and
#'         prognostication in chronic lymphocytic leukaemia. Immunogenetics. 67(1):61-6
#' }
#' 
#' @seealso Like \code{updateAlleleNames}, \link{sortAlleles} can help
#'          format a list of allele names.
#' 
#' @examples
#' # Create a vector that uses old gene/allele names.
#' alleles <- c("IGHV1-c*01", "IGHV1-f*02", "IGHV2-5*07")
#' 
#' # Update the alleles to the new names
#' updateAlleleNames(alleles)
#' 
#' @export
updateAlleleNames <- function(allele_calls) {
    . = NULL
    temporary_names = c("IGHV1-c*",
                        "IGHV1-f*",
                        "IGHV3-d*",
                        "IGHV3-h*",
                        "IGHV4-b*",
                        "IGHV5-a*",
                        "IGHV2-5*10",
                        "IGHV2-5*07")
    definitive_names = c("IGHV1-38-4*",
                         "IGHV1-69-2*",
                         "IGHV3-38-3*",
                         "IGHV3-69-1*",
                         "IGHV4-38-2*",
                         "IGHV5-10-1*",
                         "IGHV2-5*02",
                         "IGHV2-5*04")
    for (i in 1:length(temporary_names)){
        allele_calls = allele_calls %>%
            gsub(temporary_names[i], definitive_names[i], ., fixed = TRUE)
    }
    return(allele_calls)
}

#' Sort allele names
#'
#' \code{sortAlleles} returns a sorted vector of strings respresenting Ig allele
#' names. Names are first sorted by gene family, then by gene, then by allele.
#' Duplicated genes have their alleles are sorted as if they were part of their
#' non-duplicated counterparts (e.g. \code{IGHV1-69D*01} comes after \code{IGHV1-69*01} 
#' but before \code{IGHV1-69*02}), and non-localized genes (e.g. \code{IGHV1-NL1*01}) 
#' come last within their gene family.
#' 
#' @param    allele_calls  vector of strings respresenting Ig allele names.
#' @param    method        string defining the method to use when sorting alleles.
#'                         If \code{"name"} then sort in lexicographic order. If
#'                         \code{"position"} then sort by position in the locus, as
#'                         determined by the final two numbers in the gene name.
#' @return   A sorted vector of strings respresenting Ig allele names.
#' 
#' @seealso Like \code{sortAlleles}, \link{updateAlleleNames} can help
#'          format a list of allele names.
#' 
#' @examples
#' # Create a list of allele names
#' alleles <- c("IGHV1-69D*01","IGHV1-69*01","IGHV1-2*01","IGHV1-69-2*01",
#'              "IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
#'              "IGHV1-2", "IGHV1-2*02", "IGHV1-69*02")
#' 
#' # Sort the alleles by name
#' sortAlleles(alleles)
#' 
#' # Sort the alleles by position in the locus
#' sortAlleles(alleles, method="pos")
#' 
#' @export
sortAlleles <- function(allele_calls, method=c("name", "position")) { 
    # Check arguments
    method <- match.arg(method)
    
    # Standardize format of submitted alleles, first
    SUBMITTED_CALLS <- getAllele(allele_calls, first = FALSE, strip_d= FALSE)
    allele_df <- data.frame(
        list("SUBMITTED_CALLS"=SUBMITTED_CALLS,
             "SUBMITTED_NAMES"=allele_calls),
        stringsAsFactors = FALSE) %>%
        arrange(SUBMITTED_CALLS) %>%
        # Determine the family
        dplyr::mutate(FAMILY = getFamily(!! rlang::sym("SUBMITTED_CALLS"))) %>%
        # Determine the gene (exclude family); convert letters to numbers for sort
        dplyr::mutate(GENE = getGene(!! rlang::sym("SUBMITTED_CALLS"))) %>%
        dplyr::mutate(GENE1 = gsub("[^-]+[-S]([^-\\*D]+).*","\\1",!! rlang::sym("SUBMITTED_CALLS"))) %>%
        dplyr::mutate(GENE1 = as.numeric(gsub("[^0-9]+", "99", !!rlang::sym("GENE1")))) %>%
        # If there is a second gene number, determine that, too
        dplyr::mutate(GENE2 = gsub("[^-]+[-S][^-]+-?","",!! rlang::sym("GENE"))) %>%
        dplyr::mutate(GENE2 = as.numeric(gsub("[^0-9]+", "99", !!rlang::sym("GENE2")))) %>%
        dplyr::mutate(ALLELE = getAllele(!!rlang::sym("SUBMITTED_CALLS"))) %>%      
        dplyr::mutate(ALLELE = sub("[^\\*]+\\*|[^\\*]+$","", !!rlang::sym("ALLELE"))) %>%
        dplyr::mutate(ALLELE = as.numeric(sub("[[:alpha:]]+","",sub("_.+$", "", !!rlang::sym("ALLELE")))))
    
    # Convert missing values to 0, sort data frame
    allele_df[is.na(allele_df)] <- 0
    if (method == "name") {  
        sorted_df <- dplyr::arrange(allele_df, !!!rlang::syms(c("FAMILY", "GENE1", "GENE2", "ALLELE")))
    } else if (method == "position") {
        sorted_df <- dplyr::arrange(allele_df, desc(!!rlang::sym("GENE1")),
                             desc(!!rlang::sym("GENE2")), 
                             desc(!!rlang::sym("FAMILY")),
                             desc(!!rlang::sym("ALLELE")))
    }
    
    return(sorted_df$SUBMITTED_NAMES)
}

#' Clean up nucleotide sequences
#'
#' \code{cleanSeqs} capitalizes nucleotides and replaces all characters 
#' besides \code{c("A", "C", "G", "T", "-", ".")} with \code{"N"}. 
#' 
#' @param    seqs  vector of nucleotide sequences.
#' 
#' @return   A modified vector of nucleotide sequences.
#' 
#' @seealso \link{sortAlleles} and \link{updateAlleleNames} can
#'          help format a list of allele names.
#' 
#' @examples
#' # Clean messy nucleotide sequences
#' seqs <- c("AGAT.taa-GAG...ATA", "GATACAGTXXZZAGNNPPACA")
#' cleanSeqs(seqs)
#' 
#' @export
cleanSeqs <- function(seqs) {
    # . = NULL
    # seqs %>%
    #     toupper %>%
    #     gsub(".", "-", . , fixed = TRUE) %>%
    #     gsub("[^ACGT-]", "N", .) %>%
    #     return
    clean_seqs <- stri_replace_all_regex(stri_trans_toupper(seqs), "[^ACGT\\.\\-]", "N")
    if (!is.null(names(seqs))) {
        names(clean_seqs) <- names(seqs)
    }
    clean_seqs
}


# Private Functions -------------------------------------------------------

# Find muations-by-position compared to a germline
#
# \code{positionMutations} duplicates the rows of a data frame for each
# position to be analyzed and determines if each sample is mutated at that
# position
# 
# @param  data          \code{data.frame} containing repertoire data. See
#                       \link{findNovelAlleles} for a list of required
#                       columns.
# @param  germline      the germline to which all the sequences should be
#                       compared
# @param  pos_range     the range of positions within the sequence for which
#                       the rows should be duplicated and checked for mutation
# @param  seq           name of the column in \code{data} with the 
#                       aligned, IMGT-numbered, V(D)J nucleotide sequence.
#                       Default is \code{sequence_alignment} 
# @return  A data frame with rows duplicated for all the positions to be
# analyzed and a column indicating whether the position is mutated in
# comparison to the germline
#
positionMutations <- function(data, germline, pos_range, seq="sequence_alignment", v_seq_length = "v_germline_end"){
    . = NULL
    pos_db = pos_range %>%
        length() %>%
        rep("data", .) %>%
        paste(collapse=",") %>%
        paste("bind_rows(",., ")") %>%
        parse(text=.) %>%
        eval()
    pos_db$POSITION = c(sapply(pos_range, rep, nrow(data)))
    # Find which positions are mutated
    pos_db <- pos_db %>% rowwise() %>%
        mutate(NT = substring(!!rlang::sym(seq),
                              !!rlang::sym("POSITION"),
                              !!rlang::sym("POSITION"))) %>%
        mutate(GERM_NT = substring(germline, !!rlang::sym("POSITION"), !!rlang::sym("POSITION"))) %>%
        mutate(MUTATED = (!!rlang::sym("NT") != !!rlang::sym("GERM_NT") & 
                          !!rlang::sym("NT") != "N" & 
                          !!rlang::sym("NT") != "-" & 
                          !!rlang::sym("NT") != "." & 
                          !!rlang::sym("NT") != "")) %>%
        mutate(OBSERVED = (!!rlang::sym("NT") != "-" & 
                           !!rlang::sym("NT") != "." & 
                           !!rlang::sym("NT") != "")) %>% filter(!!rlang::sym("POSITION")<= !!rlang::sym(v_seq_length))
    if (any(pos_db$GERM_NT == "") && !grepl("\\.", germline) ) { 
        stop("Empty ('') GERM_NT positions found. Check you are using gapped reference germlines.")
        }
    return(pos_db)
}

# Find sequences carrying certain levels of mutation
#
# \code{mutationRangeSubset} determines the mutations in a \code{data.frame} of
# sequences and returns the subset of sequences that meet the given mutation
# count limits
# 
# @param  data          \code{data.frame} containing repertoire data. See
#                       \link{findNovelAlleles} for a list of required
#                       columns.
# @param  germline      germline to which all the sequences should be
#                       compared.
# @param  pos_range     range of positions within the sequences that should
#                       be analyzed for mutations
# @param  pos_range     range of mutation counts that sequences can have
#                       and still be included.
# @param  seq           name of the column in \code{data} with the 
#                       aligned, IMGT-numbered, V(D)J nucleotide sequence.
#                       Default is \code{sequence_alignment}.
# @return
# A data.frame containing only the subset carrying the desired levels
# of mutation
#
mutationRangeSubset <- function(data, germline, mut_range, pos_range, 
                                seq="sequence_alignment",v_seq_length = v_seq_length){
    . = NULL
    pads = paste(rep("-", min(pos_range)-1), collapse="")
    data$MUT_COUNT = data %>% rowwise() %>%
        mutate(subseq = substring(!!rlang::sym(seq), min(pos_range), min(max(pos_range),!!rlang::sym(v_seq_length)))) %>%
        pull(subseq) %>%
        paste(pads, ., sep="") %>%
        getMutatedPositions(germline) %>%
        sapply(length)
    data = data %>%
        filter(!!rlang::sym("MUT_COUNT") %in% mut_range)
    return(data)
}

# Find lower range of y-intercept confidence interval
#
# \code{findLowerY} finds the lower range of y-intercept confidence interval
# 
# @details  If mut_min is 1, a y-intercept will be searched for at 0. If
# mut_min is above 1, then the "y-intercept" will be found at x = mut_min - 1.
#
# @param    x         vector of x values.
# @param    y         vector of y values.
# @param    mut_min   value where the the lowest mutation count should be
#                     found. See details.
# @param    alpha     alpha cutoff the be used in constructing the
#                     confidence interval
#
# @return  A data frame containing only the subset carrying the desired levels
# of mutation
#
findLowerY = function(x, y, mut_min, alpha, weights){
    y = y + 1 - mut_min
    lowerY = suppressWarnings(confint(lm(x ~ y), level=1 - 2*alpha)[[1]]) #, weights = weights
    return(lowerY)
}

# Enchanced substring extraction
#
# \code{superSubstring} is an enahnced version of \code{substring} in that
# it can find disjoint positions in one call.
#
# @param    string      single string.
# @param    positions   the positions to be extracted.
#
# @return   Substring.
superSubstring = function(string, positions){
    if(length(string) != 1){ stop("Please submit only one string.") }
    chars = sapply(positions, function(x) substring(string, x, x))
    return(paste(chars, collapse=""))
}


# Layout multiple ggplots
#
# \code{multiplot} is a function provided by http://www.cookbook-r.com/ which
# allows for plotting multiple ggplot objects in one panel.
# 
# @param    ...       ggplot2 object(s).
# @param    plotlist  a list alternative to (...).
# @param    file      an unused parameter, but present in the provided function.
# @param    cols      Number of columns in layout.
# @param    layout    matrix specifying the layout. If present, \code{cols} is
#                     ignored.
# @param    heights   numeric vector A numeric vector or unit object 
#                     describing the heights of the rows in the layout. Will
#                     be passed to grid.layout. Default is all plots have 
#                     the same height.
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, heights=NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots <- length(plots)
    ncol <- cols
    nrow <- ceiling(numPlots/cols)
    if (is.null(heights)) { heights = rep(1,nrow) }
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * nrow),
                         ncol = cols, nrow = nrow)
    }
    grob <- gridExtra::arrangeGrob(grobs=plots, 
                                   nrow=nrow, ncol=ncol, layout_matrix = layout,
                                   heights=heights)
    p <- ggplot() + 
        layer(data = data.frame(x = NA),
              stat = StatIdentity, 
              position = PositionIdentity, 
              # geom = GeomDrawGrob, 
              geom = GeomCustomAnn,
              inherit.aes = FALSE, 
              params = list(grob = grob, 
                            xmin = 0,
                            xmax = 1, 
                            ymin = 0, 
                            ymax = 1)) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))
    p
}


#' Subsample repertoire
#'
#' \code{subsampleDb} will sample the same number of sequences for each gene, family
#' or allele (specified with \code{mode}) in \code{data}. Samples or subjects can
#' be subsampled indepently by setting \code{group}.
#' 
#' \code{data} will be split into gene, allele or family subsets (\code{mode}) from
#' which the same number of sequences will be subsampled. If \code{mode=gene},
#' for each gene in the field \code{gene} from \code{data}, a maximum of 
#' \code{max_n} sequences will be subsampled. Input sequences
#' that have multiple gene calls (ties), can be subsampled from any of their calls, 
#' but these duplicated samplings will be removed, and the final 
#' subsampled \code{data} will contain unique rows.
#' 
#' @param   data   \code{data.frame} containing repertoire data.
#' @param   gene   name of the column in \code{data} with allele calls. Default
#'                 is \code{v_call}.
#' @param   mode   one of \code{c("gene", "family", "allele")} defining the degree of
#'                 specificity regarding allele calls when subsetting sequences.
#'                 Determines how \code{data} will be split into subsets from 
#'                 which the same number of sequences will be subsampled. See 
#'                 also \code{group}.
#' @param   min_n  minimum number of observations to sample from each groupe. A group with 
#'                 less observations than the minimum is excluded. 
#' @param   max_n  maximum number of observations to sample for all \code{mode} groups.
#'                 If \code{NULL}, it will be set automatically to the size of 
#'                 the smallest group. If \code{max_n} is larger than the availabe 
#'                 number of sequences for any \code{mode} group, if will be 
#'                 automatically adjusted and the efective \code{max_n} used 
#'                 will be the size of the smallest \code{mode} group.
#' @param   group  columns containing additional grouping variables, e.g. sample_id.
#'                 These groups will be subsampled independently. If
#'                 \code{max_n} is \code{NULL}, a \code{max_n} will be 
#'                 automatically set for each \code{group}.
#' 
#' @return  Subsampled version of the input \code{data}.
#' 
#' @seealso \link{selectNovel}
#' 
#' @examples
#' subsampleDb(AIRRDb)
#' 
#' @export
subsampleDb <- function(data, gene="v_call", mode=c("gene", "allele", "family"), 
                        min_n=1, max_n=NULL, group=NULL) {
    
    mode <- match.arg(mode)
    
    # Check columns exist(can be NULL)
    check <- checkColumns(data, c(gene, group))
    if (check != TRUE) { stop(check) }
    
    # If group is set, 
    # split by group, then call this same function without grouping, 
    # to apply the subsample within each group independently
    if (!is.null(group)) {
        group_id <- data %>%
            group_by(!!rlang::sym(group)) %>%
            group_indices()
        ss_data <- bind_rows(lapply(split(data, group_id), subsampleDb,
                             gene=gene, mode=mode, min_n=min_n, max_n=max_n, group=NULL)
        ) 
        return(ss_data)
    }
    
    # Extract gene, allele or family assignments
    gene_func <- switch(mode,
                        allele=getAllele,
                        gene=getGene,
                        family=getFamily)
    
    # Get temporary 'gene' names (SubSet Gene) used to 
    # create subsampling groups
    data[["SS_GENE"]] <- gene_func(data[[gene]], first=FALSE)
    
        
    # Find which rows' calls contain which germline alleles
    genes <- unique(unlist(strsplit(data[["SS_GENE"]], ",")))
    allele_groups <- sapply(genes, grep, data[["SS_GENE"]],
                            fixed=TRUE, simplify=FALSE)
    allele_groups_sizes <- sapply(allele_groups, length)
    allele_groups <- allele_groups[ allele_groups_sizes >= min_n ]
    if(length(allele_groups) == 0){
        stop("Not enough sample sequences (min_n) were assigned to any gene. ",
                "Returned `NULL`")
    }
    
    # Determine the number of sequences to subsample
    if (!is.null(max_n)) {
        n <- min(allele_groups_sizes, max_n)
    } else {
        n <- min(allele_groups_sizes)
    }
    
    # Get indices subsampled sequences
    ss_idx <- lapply(allele_groups, sample, size=n, replace=FALSE, prob=NULL)
    ss_idx <- unique(unlist(ss_idx))
    
    data[ss_idx, ] %>% select(-!!rlang::sym("SS_GENE"))
}
