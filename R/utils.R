
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


#' Wrapper to check if command is found in PATH
#' 
#' @param cmd Command to test.
#' @param args Arguments for command.
#'
#' @return Logical indicating whether the command is in PATH or not.
#' @noRd
is_valid <- function(cmd = NULL, args = NULL) {
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
  valid <- is_valid(cmd = "hmmsearch", args = "-h")
  return(valid)
}


