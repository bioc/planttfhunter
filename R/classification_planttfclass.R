
#' Add Plant-TFClass classifications based on TAPscan classifications
#'
#' Plant-TFClass families, classes, and superclasses are inferred from
#' TAPscan families and subfamilies using a correspondence table between
#' the two classification systems. When a TAPscan family corresponds to
#' more than one Plant-TFClass family (e.g., TAPscan's AP2 includes
#' Plant-TFClass's AP2, ERF/DREB, and RAV), PlantTFDB subfamilies are used
#' to resolve the ambiguity. If the ambiguity cannot be resolved, all possible
#' classifications are reported, separated by ';'.
#'
#' @param families A data frame with variables \strong{Gene},
#' \strong{PlantTFDB_subfamily}, \strong{TAPscan_family}, and
#' \strong{TAPscan_subfamily}.
#'
#' @return The input data frame with three additional variables:
#' \strong{PlantTFClass_superclass}, \strong{PlantTFClass_class},
#' and \strong{PlantTFClass_family}. Classes without families have the
#' class name as family.
#' @noRd
.classify_planttfclass <- function(families = NULL) {

    map <- planttfclass_map
    map_key <- paste(map$tapscan_family, map$tapscan_subfamily, sep = "\r")
    gene_key <- paste(
        families$TAPscan_family, families$TAPscan_subfamily, sep = "\r"
    )

    # Classify each unique combination of TAPscan and PlantTFDB (sub)families
    ptfdb <- families$PlantTFDB_subfamily
    combo <- paste(gene_key, ptfdb, sep = "\t")
    first <- !duplicated(combo)
    res <- vapply(which(first), function(i) {
        .resolve_planttfclass(map[map_key == gene_key[i], ], ptfdb[i])
    }, character(3))
    res <- matrix(res, ncol = 3, byrow = TRUE)
    res <- res[match(combo, combo[first]), , drop = FALSE]

    families$PlantTFClass_superclass <- res[, 3]
    families$PlantTFClass_class <- res[, 2]
    families$PlantTFClass_family <- res[, 1]
    return(families)
}


#' Choose Plant-TFClass classification among candidates
#'
#' @param candidates A data frame with candidate Plant-TFClass
#' classifications for a gene (i.e., rows of the correspondence table between
#' TAPscan and Plant-TFClass).
#' @param planttfdb_subfamily Character scalar with the gene's PlantTFDB
#' subfamily (or subfamilies, separated by ';'), or NA.
#'
#' @return A character vector of length 3 with Plant-TFClass family,
#' class, and superclass. If more than one candidate remains after
#' resolution with PlantTFDB subfamilies, values are separated by ';' (in
#' alphabetical order of families).
#' @noRd
.resolve_planttfclass <- function(candidates, planttfdb_subfamily = NA) {

    if(nrow(candidates) > 1) {
        ptfdb <- unlist(strsplit(planttfdb_subfamily, ";"))
        keep <- !is.na(candidates$planttfdb_subfamily) &
            candidates$planttfdb_subfamily %in% ptfdb
        if(any(keep)) { candidates <- candidates[keep, ] }
    }
    candidates <- candidates[order(candidates$family, method = "radix"), ]

    collapse <- function(x) {
        x <- unique(x[!is.na(x)])
        if(length(x) == 0) { return(NA_character_) }
        return(paste(x, collapse = ";"))
    }
    res <- c(
        collapse(candidates$family), collapse(candidates$class),
        collapse(candidates$superclass)
    )
    return(res)
}
