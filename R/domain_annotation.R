

#' Annotate proteins sequences with PFAM domains
#' 
#' PFAM domains are assigned to each sequence using HMMER.
#' 
#' @param seq An AAStringSet object as returned 
#' by \code{Biostrings::readAAStringSet()}. The sequences in this object must 
#' represent only the translated sequences of primary (or longest) transcripts.
#' @param mode Character specifying how to annotate PFAM domains. One of
#' 'web' (default) or 'local'. If 'local', you must have hmmer installed in
#' your machine, and the command 'hmmsearch' must be in your PATH.
#' @param evalue Numeric indicating the E-value threshold for domain annotation
#' with hmmsearch. Only valid if parameter mode = 'local'. Default: 1e-05.
#'
#' @return A 2-column data frame with the variables \strong{Gene} 
#' and \strong{Domain}, which contain gene IDs and domain IDs, respectively.
#' @importFrom Biostrings writeXStringSet
#' @importFrom bio3d read.fasta hmmer
#' @export
#' @rdname annotate_pfam
#' @examples 
#' data(gsu)
#' seq <- gsu[1:5]
#' annotate_pfam(seq)
annotate_pfam <- function(seq = NULL, mode = "web", evalue = 1e-05) {
  seq_path <- paste0(tempdir(), "/seq.fasta")
  
  gene_names <- names(seq)
  
  if(mode == "local") {
    if(!hmmer_is_installed()) {
      stop("Could not find HMMER. Check if it is installed and in your PATH.
       Alternatively, set mode = 'web' to use the HMMER web server.")
    }
    Biostrings::writeXStringSet(seq, filepath = seq_path)
    profiles <- system.file("extdata", package = "tfhunter")
    profiles <- list.files(profiles, pattern = ".hmm", full.names = TRUE)
    annotation <- lapply(profiles, function(x) {
      out_file <- paste0("--domtblout ", tempdir(), "/search.txt")
      system2("hmmsearch", args = c(out_file, x, seq_path), stdout = FALSE)
      result <- read_hmmsearch(path = paste0(tempdir(), "/search.txt"))
      annot <- NULL
      if(!is.null(result)) {
        result <- result[result$sequence_evalue < evalue, ]
        if(nrow(result) != 0) {
          annot <- data.frame(
            Gene = result$domain_name, 
            Domain = result$query_name
          )
        }
      }
      return(annot)
    })
    annotation <- Reduce(rbind, annotation)
  } else if(mode == "web") {
    annotation <- lapply(seq_along(seq), function(x) {
      Biostrings::writeXStringSet(seq[x], filepath = seq_path)
      protein <- bio3d::read.fasta(seq_path)
      pfam <- tryCatch(
        bio3d::hmmer(protein, type="hmmscan", db="pfam", 
                     verbose = FALSE)$hit.tbl$acc,
        error = function(e) return(NA), warning = function(w) return(NA)
      )
      annot <- data.frame(
        Gene = rep(gene_names[x], length(pfam)), Domain = pfam
      )
      return(annot)
    })
    annotation <- Reduce(rbind, annotation)
  } else {
    stop("Choose one of 'web' or 'local' for parameter mode.")
  }
  return(annotation)
}

