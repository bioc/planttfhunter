
#' Annotate protein sequences with domains used for TF classification
#'
#' Protein sequences are scanned with HMMER for domains used by the
#' classification schemes of \strong{PlantTFDB} and \strong{TAPscan}. 
#' For PlantTFDB, PFAM and self-built profile HMMs are used with the 
#' domain-specific score cutoffs defined by PlantTFDB (or an E-value threshold, 
#' for domains without cutoffs). 
#' For TAPscan, TAPscan v4 profile HMMs are used with their
#' gathering (GA) thresholds, as in the original TAPscan implementation.
#'
#' @param seq An `AAStringSet` object. The sequences in this object must
#' represent only the translated sequences of primary (or longest) transcripts.
#' To create such an object from a FASTA file, 
#' use \code{Biostrings::readAAStringSet()}.
#' @param evalue Numeric indicating the E-value threshold for \code{hmmsearch}
#' to be used for PlantTFDB domains without pre-defined domain cutoffs.
#' Default: 1e-05.
#' @param threads Numeric indicating the number of threads to be used
#' by \code{hmmsearch}. Default: 1.
#'
#' @return A list of 2 data frames named \strong{PlantTFDB} and
#' \strong{TAPscan}, with domain annotation for each classification scheme.
#' Both data frames have one row per domain hit and the following variables:
#' \describe{
#'   \item{Gene}{Character, gene ID.}
#'   \item{Domain}{Character, domain ID or name.}
#'   \item{qlen}{Numeric, length of the profile HMM.}
#'   \item{c_evalue}{Numeric, conditional E-value of the domain hit.}
#'   \item{hmm_from, hmm_to}{Numeric, start and end coordinates of the alignment
#'   in the profile HMM.}
#' }
#' Rows are kept in the same order as in the HMMER output, as this order
#' is used by the TAPscan classification algorithm.
#'
#' @export
#' @rdname annotate_domains
#' @examples
#' data(gsu)
#' seq <- gsu[1:5]
#' if(hmmer_is_installed()) {
#'     annotate_domains(seq)
#' }
annotate_domains <- function(seq = NULL, evalue = 1e-05, threads = 1) {

    seq_path <- .write_seqs(seq)
    on.exit(unlink(seq_path))

    annotation <- list(
        PlantTFDB = .annotate_planttfdb(seq_path, evalue, threads),
        TAPscan = .annotate_tapscan(seq_path, threads)
    )
    return(annotation)
}


#' Write sequences to a temporary FASTA file after checking input
#'
#' @param seq An `AAStringSet` object.
#'
#' @return Path to the temporary FASTA file.
#' @importFrom Biostrings writeXStringSet
#' @importFrom methods is
#' @noRd
.write_seqs <- function(seq = NULL) {

    if(!hmmer_is_installed()) {
        stop("Could not find HMMER. Check if it is installed and in your PATH.")
    }
    if(!is(seq, "AAStringSet")) {
        stop("Argument to 'seq' must be an AAStringSet object.")
    }

    seq_path <- tempfile(pattern = "seq", fileext = ".fasta")
    Biostrings::writeXStringSet(seq, filepath = seq_path)
    return(seq_path)
}


#' Annotate sequences in a FASTA file with PlantTFDB domains
#'
#' @param seq_path Path to FASTA file.
#' @inheritParams annotate_domains
#'
#' @return A data frame as in the element \strong{PlantTFDB} of the output
#' of \code{annotate_domains()}.
#' @noRd
.annotate_planttfdb <- function(seq_path, evalue = 1e-05, threads = 1) {

    hmms <- file.path(
        system.file("extdata", package = "planttfhunter"),
        c("PFAM.hmm", "self_built.hmm")
    )
    hits <- lapply(hmms, function(x) {
        hits <- .run_hmmsearch(seq_path, x, threads = threads)
        .filter_hmmer(hits, evalue)
    })
    hits <- do.call(rbind, hits)

    return(.format_hits(hits))
}


#' Annotate sequences in a FASTA file with TAPscan domains
#'
#' @param seq_path Path to FASTA file.
#' @inheritParams annotate_domains
#'
#' @return A data frame as in the element \strong{TAPscan} of the output
#' of \code{annotate_domains()}.
#' @noRd
.annotate_tapscan <- function(seq_path, threads = 1) {

    hits <- .run_hmmsearch(
        seq_path, .tapscan_hmm_path(), extra_args = "--cut_ga",
        threads = threads
    )

    return(.format_hits(hits))
}


#' Format hmmsearch hits as domain annotation
#'
#' @param hits A data frame with hmmsearch hits as returned
#' by \code{.read_hmmsearch()}.
#'
#' @return A data frame with variables \strong{Gene}, \strong{Domain},
#' \strong{qlen}, \strong{c_evalue}, \strong{hmm_from}, and
#' \strong{hmm_to}, as in the output of \code{annotate_domains()}.
#' @noRd
.format_hits <- function(hits = NULL) {

    annotation <- data.frame(
        Gene = hits$target_name,
        Domain = hits$query_name,
        qlen = hits$qlen,
        c_evalue = hits$c_evalue,
        hmm_from = hits$hmm_from,
        hmm_to = hits$hmm_to
    )
    return(annotation)
}
