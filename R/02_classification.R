
#' Identify TFs and classify them in families
#'
#' TFs are classified using the classification schemes of PlantTFDB,
#' TAPscan, and Plant-TFClass.
#'
#' @param domain_annotation A list of data frames with domain annotation,
#' as returned by \code{annotate_domains()}.
#' @param scheme Character indicating which classification scheme(s) to
#' return. One of 'all' (all schemes), 'PlantTFDB', 'TAPscan', or
#' 'PlantTFClass'. Default: 'all'.
#'
#' @return A data frame with one row per gene and the following variables
#' (if \strong{scheme = "all"}):
#' \describe{
#'   \item{Gene}{Gene ID.}
#'   \item{PlantTFDB_family}{TF family according to PlantTFDB.}
#'   \item{PlantTFDB_subfamily}{TF subfamily according to PlantTFDB.}
#'   \item{TAPscan_family}{TF family according to TAPscan.}
#'   \item{TAPscan_subfamily}{TF subfamily according to TAPscan.}
#'   \item{TAP_class}{Class of transcription-associated protein (TAP)
#'   according to TAPscan: transcription factor (TF), transcriptional
#'   regulator (TR), or putative TAP (PT). See \code{data(tap_class)}.}
#'   \item{PlantTFClass_superclass}{TF superclass according to
#'   Plant-TFClass.}
#'   \item{PlantTFClass_class}{TF class according to Plant-TFClass.}
#'   \item{PlantTFClass_family}{TF family according to Plant-TFClass.}
#' }
#' Only genes classified by at least one scheme are included. Genes not
#' classified by a scheme have NA in the respective variables. Genes
#' assigned to more than one PlantTFDB family have families (and
#' subfamilies, in the same order) separated by ';'.
#'
#' If \strong{scheme} is 'PlantTFDB', 'TAPscan', or 'PlantTFClass', only
#' variable \strong{Gene} and the variables of the chosen scheme are
#' returned, and only genes classified by that scheme are included. Note that
#' TAPscan classifies not only TFs, but also transcriptional regulators and
#' putative TAPs, which can be distinguished using the variable
#' \strong{TAP_class}.
#' 
#' @details
#' PlantTFDB and TAPscan classifications are based on the domains found in
#' each gene. Plant-TFClass is a structure-based classification of TFs
#' into superclasses, classes, and families, and Plant-TFClass
#' classifications are inferred from TAPscan classifications using a
#' correspondence table between the two classification systems. When a
#' TAPscan subfamily corresponds to more than one Plant-TFClass family
#' (e.g., TAPscan's AP2 includes Plant-TFClass's AP2, ERF/DREB, and RAV),
#' PlantTFDB subfamilies are used to resolve the ambiguity. If the ambiguity
#' cannot be resolved (e.g., the gene was not classified by PlantTFDB),
#' all possible classifications are reported, separated by ';'.
#'
#' In all schemes, classification levels without a more specific level
#' are filled with the name of the more general level. For example,
#' PlantTFDB family bHLH has no subfamilies, so its subfamily is bHLH.
#'
#' The classification schemes are available in the data sets
#' \code{planttfdb_scheme} (PlantTFDB), \code{tapscan_scheme}
#' (TAPscan), and \code{planttfclass_scheme} (correspondence between
#' TAPscan and Plant-TFClass).
#'
#' @references
#' Jin, J., Tian, F., Yang, D. C., Meng, Y. Q., Kong, L., Luo, J., &
#' Gao, G. (2017). PlantTFDB 4.0: toward a central hub for transcription
#' factors and regulatory interactions in plants. Nucleic Acids Research,
#' 45(D1), D1040-D1045.
#'
#' Petroll, R., Varshney, D., Hiltemann, S., Finke, H., Schreiber, M.,
#' de Vries, J., & Rensing, S. A. (2025). Enhanced sensitivity of TAPscan v4
#' enables comprehensive analysis of streptophyte transcription factor
#' evolution. The Plant Journal, 121(1), e17184.
#'
#' Blanc-Mathieu, R., Dumas, R., Turchi, L., Lucas, J., & Parcy, F. (2024).
#' Plant-TFClass: a structural classification for plant transcription
#' factors. Trends in Plant Science, 29(1), 40-51.
#'
#' @export
#' @rdname classify_tfs
#' @examples
#' data(gsu_annotation)
#' families <- classify_tfs(gsu_annotation)
#' head(families)
#'
#' # Only TAPscan classifications, with TAP classes (TF, TR, or PT)
#' tapscan_families <- classify_tfs(gsu_annotation, scheme = "TAPscan")
#' table(tapscan_families$TAP_class)
classify_tfs <- function(
        domain_annotation = NULL,
        scheme = c("all", "PlantTFDB", "TAPscan", "PlantTFClass")
) {

    scheme <- match.arg(scheme)
    .check_domain_annotation(domain_annotation)

    ptfdb <- .classify_planttfdb(domain_annotation$PlantTFDB)
    names(ptfdb) <- c("Gene", .scheme_levels("PlantTFDB"))
    tapscan <- .classify_tapscan(domain_annotation$TAPscan)
    names(tapscan) <- c("Gene", .scheme_levels("TAPscan"))

    families <- merge(ptfdb, tapscan, by = "Gene", all = TRUE)
    families <- families[order(families$Gene, method = "radix"), ]
    families <- .add_tap_class(families)
    families <- .classify_planttfclass(families)

    # Select variables and genes of the chosen scheme(s)
    if(scheme == "all") {
        cols <- unlist(lapply(
            c("PlantTFDB", "TAPscan", "PlantTFClass"), .scheme_columns
        ))
    } else {
        cols <- .scheme_columns(scheme)
        families <- families[!is.na(families[[cols[1]]]), ]
    }
    families <- families[, c("Gene", cols)]
    rownames(families) <- NULL
    
    return(families)
}


#' Check domain annotation passed to \code{classify_tfs()}
#'
#' @inheritParams classify_tfs
#'
#' @return TRUE (invisibly) if the input is valid. Otherwise, an error is
#' thrown.
#' @importFrom methods is
#' @noRd
.check_domain_annotation <- function(domain_annotation) {

    cols <- c("Gene", "Domain", "qlen", "c_evalue", "hmm_from", "hmm_to")
    valid <- is(domain_annotation, "list") &&
        all(c("PlantTFDB", "TAPscan") %in% names(domain_annotation)) &&
        all(vapply(domain_annotation[c("PlantTFDB", "TAPscan")], function(x) {
            is(x, "data.frame") && all(cols %in% names(x))
        }, logical(1)))

    if(!valid) {
        stop(
            "Argument to 'domain_annotation' must be a list of data ",
            "frames named 'PlantTFDB' and 'TAPscan', as returned by ",
            "annotate_domains()."
        )
    }
    return(invisible(TRUE))
}


#' Get variables with classification levels of a classification scheme
#'
#' @param scheme Character indicating the classification scheme. One of
#' 'PlantTFDB', 'TAPscan', or 'PlantTFClass'.
#'
#' @return A named character vector with the names of the variables in the
#' output of \code{classify_tfs()} for each classification level, from
#' the most general to the most specific level. Vector names are the
#' classification levels.
#' @noRd
.scheme_levels <- function(scheme) {

    levels <- switch(
        scheme,
        PlantTFDB = c(
            family = "PlantTFDB_family", subfamily = "PlantTFDB_subfamily"
        ),
        TAPscan = c(
            family = "TAPscan_family", subfamily = "TAPscan_subfamily"
        ),
        PlantTFClass = c(
            superclass = "PlantTFClass_superclass",
            class = "PlantTFClass_class", family = "PlantTFClass_family"
        )
    )
    return(levels)
}


#' Get all variables of a classification scheme
#'
#' @param scheme Character indicating the classification scheme. One of
#' 'PlantTFDB', 'TAPscan', or 'PlantTFClass'.
#'
#' @return A character vector with the names of the variables of the
#' classification scheme in the output of \code{classify_tfs()}, including
#' classification levels and additional variables (i.e., TAP class for
#' TAPscan).
#' @noRd
.scheme_columns <- function(scheme) {

    cols <- unname(.scheme_levels(scheme))
    if(scheme == "TAPscan") { cols <- c(cols, "TAP_class") }
    return(cols)
}
