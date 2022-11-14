

#' Read and parse hmmsearch output
#'
#' @param path Path to hmmsearch output in --domtblout format.
#'
#' @return A parsed data frame with the 23 columns HMMER outputs by default.
#' 
#' @noRd
#' @importFrom utils read.csv
read_hmmsearch <- function(path = NULL) {
  
  stab <- gsub(pattern = "  *", replacement = "\t", readLines(path))
  stab <- stab[!startsWith(stab, "#")]
  if(length(stab) != 0) {
    stab <- paste0(stab, collapse = "\n")
    stab <- utils::read.csv(textConnection(stab), sep = "\t", header = FALSE)
    colnames(stab) <- c(
      "domain_name", "domain_accession", "domain_len",
      "query_name", "query_accession", "qlen",
      "sequence_evalue", "sequence_score", "sequence_bias",
      "domain_N", "domain_of", "domain_cevalue", "domain_ievalue",
      "domain_score", "domain_bias",
      "hmm_from", "hmm_to",
      "ali_from", "ali_to",
      "env_from", "env_to",
      "acc", "description"
    )
  } else {
    stab <- NULL
  }
  return(stab)
}

#' Filter HMMER results based on domain cutoffs for each domain
#'
#' @param hmmer_results A data frame with HMMER results as returned 
#' by \code{read_hmmsearch}.
#' @param evalue Numeric indicating the E-value threshold for hmmsearch 
#' to be used for domains without pre-defined domain cutoffs. 
#' Only valid if parameter mode = 'local'. Default: 1e-05.
#' 
#' @return The same data frame passed as input, but filtered based
#' on domain cutoffs.
#' @noRd
filter_hmmer <- function(hmmer_results = NULL, evalue = 1e-05) {
    
    cutoffs <- score_cutoff
    
    res_cutoff <- merge(hmmer_results, cutoffs,
                        by.x = "query_name", by.y = "domain")
    
    # 1) Filter by domain cutoff
    res_filt1 <- res_cutoff[!is.na(res_cutoff$domaincutoff), ]
    res_filt1 <- res_filt1[res_filt1$domain_score >= res_filt1$domaincutoff, ]
    
    # 2) Filter by evalue
    res_filt2 <- res_cutoff[is.na(res_cutoff$domaincutoff), ]
    res_filt2 <- res_filt2[res_filt2$sequence_evalue < evalue, ]
    
    # Combine results
    res_filtered <- rbind(res_filt1, res_filt2) 
    return(res_filtered)
}


#' Wrapper to check if command is found in PATH
#' 
#' @param cmd Command to test.
#' @param args Arguments for command.
#'
#' @return Logical indicating whether the command is in PATH or not.
#' @noRd
is_valid_cmd <- function(cmd = NULL, args = NULL) {
  found <- tryCatch(
    system2(cmd, args = args, stdout = FALSE, stderr = FALSE), 
    error = function(e) return(FALSE),
    warning = function(w) return(FALSE)
  )
  if(!isFALSE(found)) {
    found <- TRUE
  }
  return(found)
}

#' Check if HMMER is installed
#' 
#' @return Logical indicating whether HMMER is installed or not.
#' @export
#' @rdname hmmer_is_installed
#' @examples 
#' hmmer_is_installed()
hmmer_is_installed <- function() {
  valid <- is_valid_cmd(cmd = "hmmsearch", args = "-h")
  return(valid)
}


