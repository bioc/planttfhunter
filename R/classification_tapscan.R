
#' Identify TFs and classify them in families using TAPscan's scheme
#'
#' @param domain_annotation A data frame with domain annotation, as
#' in the element \strong{TAPscan} of the output of \code{annotate_domains()}.
#'
#' @return A 3-column data frame with variables \strong{Gene},
#' \strong{Family}, and \strong{Subfamily}. Families without subfamilies
#' have the family name as subfamily.
#' @noRd
.classify_tapscan <- function(domain_annotation = NULL) {

    empty <- data.frame(
        Gene = character(0), Family = character(0), Subfamily = character(0)
    )
    hits <- .tapscan_filter_coverage(domain_annotation)
    if(nrow(hits) == 0) { return(empty) }

    hits <- .tapscan_collapse_hits(hits)
    hits <- .tapscan_exclude_similar(hits)
    tokens <- .tapscan_tokens(hits)

    fams <- .tapscan_assign_families(tokens)
    if(nrow(fams) == 0) { return(empty) }

    fams <- cbind(Gene = fams$Gene, .tapscan_family_map(fams$Family))
    rownames(fams) <- NULL
    return(fams)
}


#' Filter TAPscan domain hits by domain-specific coverage cutoffs
#'
#' @param domain_annotation A data frame with domain annotation, as
#' in the element \strong{TAPscan} of the output of \code{annotate_domains()}.
#'
#' @return A data frame with variables \strong{Gene}, \strong{Domain}, and
#' \strong{c_evalue} for hits whose coverage of the profile HMM is greater
#' than the domain-specific cutoff.
#' @noRd
.tapscan_filter_coverage <- function(domain_annotation = NULL) {

    hits <- domain_annotation
    coverage <- (hits$hmm_to - hits$hmm_from + 1) / hits$qlen
    cutoff <- tapscan_coverage$cutoff[
        match(hits$Domain, tapscan_coverage$domain)
    ]
    cutoff[is.na(cutoff)] <- 0

    hits <- hits[coverage > cutoff, c("Gene", "Domain", "c_evalue")]
    return(hits)
}


#' Collapse consecutive hits of the same domain in the same gene
#'
#' @param hits A data frame as returned by \code{.tapscan_filter_coverage()}.
#'
#' @return A data frame with variables \strong{Gene}, \strong{Domain},
#' \strong{c_evalue}, and \strong{count} (number of hits), with one row per
#' gene-domain pair, sorted by gene and E-value.
#' @noRd
.tapscan_collapse_hits <- function(hits = NULL) {

    n <- nrow(hits)
    key <- paste0(hits$Gene, hits$Domain)
    run_start <- c(TRUE, key[-1] != key[-n])
    run_end <- which(c(run_start[-1], TRUE))
    count <- diff(c(0L, run_end))

    ev <- hits$c_evalue
    prev_end <- pmax(run_end - 1L, 1L)
    run_ev <- ifelse(count >= 2, pmax(ev[run_end], ev[prev_end]), ev[run_end])

    collapsed <- data.frame(
        Gene = hits$Gene[run_end],
        Domain = hits$Domain[run_end],
        c_evalue = run_ev,
        count = count
    )
    collapsed <- collapsed[order(
        collapsed$Gene, collapsed$c_evalue, method = "radix"
    ), ]
    rownames(collapsed) <- NULL
    return(collapsed)
}


#' Remove hits to domains that are similar to better-scored domains
#'
#' For groups of similar domains (e.g., Myb_DNA-binding and G2-like_Domain),
#' only the first domain found for each gene (i.e., the one with the lowest
#' E-value) is kept. Domain FIE_clipped_for_HMM is removed if WD40 was found
#' before it.
#'
#' @param hits A data frame as returned by \code{.tapscan_collapse_hits()}.
#'
#' @return The input data frame without hits to similar domains.
#' @importFrom stats ave
#' @noRd
.tapscan_exclude_similar <- function(hits = NULL) {

    line <- paste(hits$Gene, hits$Domain, hits$c_evalue, hits$count, sep = "\t")
    prev_gene <- c("", hits$Gene[-nrow(hits)])
    idx <- split(seq_along(line), prev_gene)
    same <- logical(length(line))
    same[unlist(idx)] <- unlist(lapply(seq_along(idx), function(k) {
        .regex_starts_with(line[idx[[k]]], names(idx)[k])
    }))
    block <- cumsum(!same)

    # Groups of similar domains: only the first domain of each group is kept
    groups <- c(
        `Myb_DNA-binding` = "myb", `G2-like_Domain` = "myb",
        PHD = "phd", `Alfin-like` = "phd", C1_2 = "phd",
        GATA = "gata", `zf-Dof` = "gata"
    )
    group <- groups[hits$Domain]
    in_group <- !is.na(group)
    gkey <- paste(block, group)[in_group]
    first_dom <- hits$Domain[in_group][match(gkey, gkey)]
    drop_similar <- logical(nrow(hits))
    drop_similar[in_group] <- hits$Domain[in_group] != first_dom

    # Remove FIE if WD40 was found before it
    wd40_before <- ave(as.integer(hits$Domain == "WD40"), block, FUN = cumsum)
    drop_fie <- hits$Domain == "FIE_clipped_for_HMM" & wd40_before > 0

    hits <- hits[!drop_similar & !drop_fie, ]
    rownames(hits) <- NULL
    return(hits)
}


#' Check if strings start with a regular expression
#'
#' @param x Character vector.
#' @param pattern Character scalar with a regular expression.
#'
#' @return A logical vector indicating whether each element of \strong{x}
#' starts with \strong{pattern}. If \strong{pattern} is not a valid regular
#' expression, it is matched literally.
#' @noRd
.regex_starts_with <- function(x, pattern) {
    res <- tryCatch(
        grepl(paste0("^", pattern), x, perl = TRUE),
        error = function(e) startsWith(x, pattern),
        warning = function(w) startsWith(x, pattern)
    )
    return(res)
}


#' Get domain tokens used for TAPscan classification
#'
#' Domains are converted to tokens, which are used by TAPscan rules.
#' Tokens are the domain names themselves, plus pseudo-domains MYB-2R,
#' MYB-3R, and MYB-4R (for genes with 2, 3, or 4 Myb_DNA-binding domains),
#' and two_or_more_LIM (for genes with more than one LIM domain).
#'
#' @param hits A data frame as returned by \code{.tapscan_exclude_similar()}.
#'
#' @return A 2-column data frame with variables \strong{Gene} and
#' \strong{Token}, sorted by gene and E-value.
#' @noRd
.tapscan_tokens <- function(hits = NULL) {

    is_myb <- hits$Domain == "Myb_DNA-binding" & hits$count %in% c(2, 3, 4)
    is_lim <- hits$Domain == "LIM" & hits$count > 1
    pseudo <- rep(NA_character_, nrow(hits))
    pseudo[is_myb] <- paste0("MYB-", hits$count[is_myb], "R")
    pseudo[is_lim] <- "two_or_more_LIM"

    # Interleave pseudo-domains (when present) before their domains
    pos <- seq_len(nrow(hits)) * 2
    tokens <- data.frame(
        Gene = c(hits$Gene, hits$Gene),
        Token = c(pseudo, hits$Domain),
        pos = c(pos - 1, pos)
    )
    tokens <- tokens[!is.na(tokens$Token), ]
    tokens <- tokens[order(tokens$pos), c("Gene", "Token")]
    rownames(tokens) <- NULL
    return(tokens)
}


#' Assign genes to TAPscan families based on domain tokens
#'
#' @param tokens A data frame as returned by \code{.tapscan_tokens()}.
#'
#' @return A 2-column data frame with variables \strong{Gene} and
#' \strong{Family}, with TAPscan families as in the classification rules.
#' Only genes assigned to a family are included.
#' @noRd
.tapscan_assign_families <- function(tokens = NULL) {

    rules <- tapscan_rules
    families <- unique(rules$family)
    universe <- unique(c(tokens$Token, rules$domain))

    # Gene x token presence matrix
    genes <- unique(tokens$Gene)
    p <- unclass(table(
        factor(tokens$Gene, levels = genes),
        factor(tokens$Token, levels = universe)
    )) > 0

    # Token x family incidence matrices for 'should' and 'should not' rules
    incidence <- function(type) {
        r <- rules[rules$type == type, ]
        unclass(table(
            factor(r$domain, levels = universe),
            factor(r$family, levels = families)
        ))
    }
    s <- incidence("should")
    sn <- incidence("should not")
    n_should <- colSums(s)

    candidate <- t(t(p %*% s) == n_should) & (p %*% sn) == 0
    candidate[, n_should == 0] <- FALSE

    # Resolve genes with more than one candidate family
    should_str <- vapply(families, function(f) {
        d <- rules$domain[rules$family == f & rules$type == "should"]
        paste0(";", paste(d, collapse = ";"), ";")
    }, character(1))
    gene_tokens <- split(tokens$Token, factor(tokens$Gene, levels = genes))
    fam <- vapply(seq_along(genes), function(i) {
        cand <- families[candidate[i, ]]
        if(length(cand) == 0) { return(NA_character_) }
        .tapscan_tiebreak(cand, should_str[cand], gene_tokens[[i]])
    }, character(1))

    fam_df <- data.frame(Gene = genes, Family = fam)
    fam_df <- fam_df[!is.na(fam_df$Family), ]
    return(fam_df)
}


#' Choose a single TAPscan family for genes with multiple candidate families
#'
#' This function reproduces the tie-breaking procedure of the original code.
#' The first candidate family is chosen by default. Then, each of
#' the other candidates is compared to the first one by looking at the
#' gene's domain tokens in order. For the first token that is found in only
#' one of the two families' 'should' domains, the gene is assigned to that
#' family. If the token is found only in the candidate, the candidate
#' replaces the current choice.
#'
#' @param cand Character vector of candidate families, in order.
#' @param should_str Character vector with each candidate's 'should' domains,
#' in the format ';domain1;domain2;'.
#' @param tokens Character vector with the gene's domain tokens, in order.
#'
#' @return Character scalar with the chosen family.
#' @noRd
.tapscan_tiebreak <- function(cand, should_str, tokens) {

    if(length(cand) == 1) { return(cand) }

    in_first <- vapply(
        tokens, grepl, logical(1), x = should_str[1], fixed = TRUE
    )
    switch_to <- vapply(seq_along(cand)[-1], function(k) {
        in_cand <- vapply(
            tokens, grepl, logical(1), x = should_str[k], fixed = TRUE
        )
        decisive <- which(xor(in_first, in_cand))
        length(decisive) > 0 && in_cand[decisive[1]]
    }, logical(1))

    chosen <- if(any(switch_to)) cand[-1][max(which(switch_to))] else cand[1]
    return(chosen)
}


#' Map TAPscan families to final family and subfamily names
#'
#' @param family Character vector of families as in the classification rules.
#'
#' @return A 2-column data frame with variables \strong{Family} and
#' \strong{Subfamily}. Families without subfamilies have the family name
#' as subfamily (as in TAPscan's family statistics).
#' @noRd
.tapscan_family_map <- function(family = NULL) {

    # Families that are merged into a single family
    merged <- c(
        `GARP_ARR-B_Myb` = "GARP_ARR-B", `GARP_ARR-B_G2` = "GARP_ARR-B",
        bZIP1 = "bZIP", bZIP2 = "bZIP", bZIPAUREO = "bZIP", bZIPCDD = "bZIP",
        HRT = "ET", GIY_YIG = "ET"
    )
    # Subfamilies and the families they belong to
    parent <- c(
        C2H2 = "C2H2", C2H2_IDD = "C2H2", bHLH = "bHLH", bHLH_TCP = "bHLH",
        `NF-YA` = "NFY", `NF-YB` = "NFY", `NF-YC` = "NFY",
        C1HDZ = "HDZ", C2HDZ = "HDZ", C3HDZ = "HDZ", C4HDZ = "HDZ",
        MYST = "HAT", CBP = "HAT", TAFII250 = "HAT", GNAT = "HAT",
        LOB1 = "LBD", LOB2 = "LBD", `MYB-related` = "MYB",
        `SWI/SNF_SWI3` = "MYB-related", `MYB-2R` = "MYB", `MYB-3R` = "MYB",
        `MYB-4R` = "MYB", RKD = "RWP-RK", NLP = "RWP-RK",
        AP2 = "AP2", CRF = "AP2"
    )

    fam <- ifelse(family %in% names(merged), merged[family], family)
    has_parent <- fam %in% names(parent)
    fam_map <- data.frame(
        Family = ifelse(has_parent, parent[fam], fam),
        Subfamily = fam
    )
    return(fam_map)
}


#' Add TAP classes (TF, TR, or PT) to TAPscan classifications
#'
#' @param families A data frame with variable \strong{TAPscan_subfamily}.
#'
#' @return The input data frame with an additional variable
#' \strong{TAP_class} indicating whether each TAPscan subfamily contains
#' transcription factors (TF), transcriptional regulators (TR), or putative
#' TAPs (PT). Subfamilies without a TAP class in TAPscan's website have NA.
#' @noRd
.add_tap_class <- function(families = NULL) {

    families$TAP_class <- tap_class_map$TAP_class[
        match(families$TAPscan_subfamily, tap_class_map$TAP_family)
    ]
    return(families)
}
