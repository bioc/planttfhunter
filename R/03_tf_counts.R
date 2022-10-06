
#' Get TF frequencies for each species as a SummarizedExperiment object
#'
#' This function identifies and classifies TFs, and returns TF counts for 
#' each family as a SummarizedExperiment object
#'
#' @param proteomes List of \strong{AAStringSet} objects
#' @param species_metadata (Optional) A data frame containing species names 
#' in row names (names must match element names in 
#' the \strong{proteomes} list), and species metadata 
#' (e.g., taxonomic information, ecological information)
#' in columns. If NULL, the colData of the SummarizedExperiment object
#' will be empty.
#'
#' @return A SummarizedExperiment object containing transcription factor
#' frequencies per family in each species, as well as species metadata (if
#' \strong{species_metadata} is not NULL).
#' 
#' @export
#' @rdname get_tf_counts
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods is
#' @examples 
#' data(gsu)
#' 
#' set.seed(123)
#' # Pick random subsets of 100 genes to simulate other species
#' proteomes <- list(
#'     Gsu1 = gsu[sample(names(gsu), 50, replace = FALSE)],
#'     Gsu2 = gsu[sample(names(gsu), 50, replace = FALSE)],
#'     Gsu3 = gsu[sample(names(gsu), 50, replace = FALSE)],
#'     Gsu4 = gsu[sample(names(gsu), 50, replace = FALSE)]
#' )
#' 
#' # Create species metadata
#' species_metadata <- data.frame(
#'     row.names = names(proteomes),
#'     Division = "Rhodophyta",
#'     Origin = c("US", "Belgium", "China", "Brazil")
#' )
#' 
#' # Get SummarizedExperiment object
#' se <- get_tf_counts(proteomes, species_metadata)
get_tf_counts <- function(proteomes, species_metadata = NULL) {
    
    # Check 1: is `proteomes` a list of AAStringSet objects?
    if(!is(proteomes, "list") || !is(proteomes[[1]], "AAStringSet")) {
        stop("Argument to 'proteomes' must be a list of AAStringSet objects.")
    }
    # Check 2: handle colData depending on `species_metadata`'s class
    coldata <- data.frame(row.names = names(proteomes))
    if(is(species_metadata, "data.frame")) {
        if(!identical(sort(rownames(species_metadata)), names(proteomes))) {
            stop("Names of 'proteomes' must match rownames of 'species_metadata'.")
        }
        coldata <- species_metadata
    }
    
    # Get a list of TFs frequencies per family
    freqs <- lapply(seq_along(proteomes), function(x) {
        
        species <- names(proteomes)[x]
        domains <- annotate_pfam(proteomes[[x]])
        families <- classify_tfs(domains)
        counts <- as.data.frame(table(families$Family))
        names(counts) <- c("Family", species)
        return(counts)
    })
    # Combine list into a quantitative matrix
    freqs_mat <- Reduce(function(x, y) merge(x, y, all = TRUE), freqs)
    rownames(freqs_mat) <- freqs_mat$Family
    freqs_mat$Family <- NULL
    freqs_mat <- as.matrix(freqs_mat)
    freqs_mat[is.na(freqs_mat)] <- 0

    # Create SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = freqs_mat), colData = coldata
    )
    return(se)
}
