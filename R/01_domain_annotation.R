

#' Annotate proteins sequences with PFAM domains
#' 
#' PFAM domains are assigned to each sequence using HMMER.
#' 
#' @param seq An AAStringSet object as returned 
#' by \code{Biostrings::readAAStringSet()}. The sequences in this object must 
#' represent only the translated sequences of primary (or longest) transcripts.
#' @param evalue Numeric indicating the E-value threshold for hmmsearch 
#' to be used for domains without pre-defined domain cutoffs. 
#' Only valid if parameter mode = 'local'. Default: 1e-05.
#'
#' @return A 2-column data frame with the variables \strong{Gene} 
#' and \strong{Domain}, which contain gene IDs and domain IDs, respectively.
#' 
#' @importFrom Biostrings writeXStringSet
#' @export
#' @rdname annotate_pfam
#' @examples 
#' data(gsu)
#' seq <- gsu[1:5]
#' if(hmmer_is_installed()) {
#'     annotate_pfam(seq)
#' }
annotate_pfam <- function(seq = NULL, evalue = 1e-05) {
    
    if(!hmmer_is_installed()) {
        stop("Could not find HMMER. Check if it is installed and in your PATH.")
    }
    
    # Write file to temporary directory
    seq_path <- tempfile("seq", fileext = ".fasta")
    Biostrings::writeXStringSet(seq, filepath = seq_path)
    
    # Run hmmsearch using pre-built HMM profiles in extdata/
    profiles <- system.file("extdata", package = "planttfhunter")
    profiles <- list.files(profiles, pattern = ".hmm", full.names = TRUE)
    
    annotation <- lapply(profiles, function(x) {
        out_file <- tempfile(pattern = "search", fileext = ".txt")
        cmd_args <- c("--domtblout", out_file, x, seq_path)
        system2("hmmsearch", args = cmd_args, stdout = FALSE)
        result <- read_hmmsearch(path = out_file)
        annot <- NULL
        if(!is.null(result)) {
            result <- filter_hmmer(result)
            
            if(nrow(result) != 0) {
                annot <- data.frame(
                    Gene = result$domain_name, 
                    Domain = result$query_name
                )
            }
        }
        unlink(out_file)
        return(annot)
    })
    annotation <- Reduce(rbind, annotation)
    unlink(seq_path)
    return(annotation)
}

