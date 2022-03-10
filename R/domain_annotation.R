

#' Annotate proteins sequences with PFAM domains
#' 
#' PFAM domains are assigned to each sequence using the HMMER web server.
#' 
#' @param seq An AAStringSet object as returned 
#' by \code{Biostrings::readAAStringSet()}.
#' 
#'
#' @importFrom Biostrings writeXStringSet
#' @importFrom bio3d read.fasta hmmer
#' @export
#' @rdname annotate_pfam
#' @examples 
#' data(gsu)
#' seq <- gsu[1:5]
#' annotate_pfam(seq)
annotate_pfam <- function(seq = NULL) {
  seq_path <- paste0(tempdir(), "/seq.fasta")

  gene_names <- names(seq)
  annotation <- lapply(seq_along(seq), function(x) {
    Biostrings::writeXStringSet(seq[x], filepath = seq_path)
    protein <- bio3d::read.fasta(seq_path)
    pfam <- tryCatch(
      bio3d::hmmer(protein, type="hmmscan", db="pfam", 
                   verbose = FALSE)$hit.tbl$acc,
      error = function(e) return(NA),
      warning = function(w) return(NA)
    )
    
    annot <- data.frame(
      Gene = rep(gene_names[x], length(pfam)),
      PFAM = pfam
    )
    return(annot)
  })
  
  annotation <- Reduce(rbind, annotation)
  return(annotation)
}

