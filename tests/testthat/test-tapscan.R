#----Load data------------------------------------------------------------------
## Output of hmmsearch against TAPscan domains for proteins of Ectocarpus sp. 7
domtbl <- .read_hmmsearch(test_path("fixtures", "ectocarpus.domtblout.gz"))
ecto_annot <- data.frame(
    Gene = domtbl$target_name, Domain = domtbl$query_name, qlen = domtbl$qlen,
    c_evalue = domtbl$c_evalue, hmm_from = domtbl$hmm_from,
    hmm_to = domtbl$hmm_to
)

## Families obtained with the original TAPscan code (v4.76)
expected <- read.table(
    test_path("fixtures", "ectocarpus_tapscan_families.txt.gz"),
    sep = ";", col.names = c("Gene", "Family", "Subfamily")
)

# Families without subfamilies ('-' in TAPscan's output)
no_sub <- expected$Subfamily == "-"
expected$Subfamily[no_sub] <- expected$Family[no_sub]

#----Start tests----------------------------------------------------------------
test_that(".classify_tapscan() reproduces the original TAPscan results", {
    fams <- .classify_tapscan(ecto_annot)
    fams <- fams[order(fams$Gene, method = "radix"), ]
    expected <- expected[order(expected$Gene, method = "radix"), ]
    rownames(fams) <- rownames(expected) <- NULL

    expect_equal(fams, expected)
})


test_that(".tapscan_filter_coverage() filters hits by coverage", {
    annot <- data.frame(
        Gene = c("g1", "g2", "g3"), Domain = c("AP2", "AP2", "unknown"),
        qlen = 100, c_evalue = 1e-10, hmm_from = 1, hmm_to = c(10, 90, 1)
    )
    hits <- .tapscan_filter_coverage(annot)
    expect_equal(hits$Gene, c("g2", "g3"))
    expect_equal(names(hits), c("Gene", "Domain", "c_evalue"))
})


test_that(".tapscan_collapse_hits() collapses consecutive hits", {
    hits <- data.frame(
        Gene = c("g2", "g2", "g2", "g1", "g1"),
        Domain = c("Myb_DNA-binding", "Myb_DNA-binding", "Myb_DNA-binding",
                   "AP2", "B3"),
        c_evalue = c(1e-05, 1e-10, 1e-08, 1e-03, 1e-20)
    )
    collapsed <- .tapscan_collapse_hits(hits)
    expect_equal(collapsed$Gene, c("g1", "g1", "g2"))
    expect_equal(collapsed$Domain, c("B3", "AP2", "Myb_DNA-binding"))
    expect_equal(collapsed$count, c(1, 1, 3))
    expect_equal(collapsed$c_evalue, c(1e-20, 1e-03, 1e-08))
})


test_that(".tapscan_exclude_similar() keeps only the best similar domain", {
    hits <- data.frame(
        Gene = c("g1", "g1", "g1", "g1", "g2", "g2"),
        Domain = c("G2-like_Domain", "Myb_DNA-binding", "WD40",
                   "FIE_clipped_for_HMM", "zf-Dof", "GATA"),
        c_evalue = 1e-10, count = 1
    )
    filtered <- .tapscan_exclude_similar(hits)
    expect_equal(filtered$Domain, c("G2-like_Domain", "WD40", "zf-Dof"))

    # Flags are not reset if gene ID starts with the previous gene ID
    hits$Gene <- c("g1", "g1", "g1", "g1", "g10", "g10")
    hits$Domain[5:6] <- c("Myb_DNA-binding", "PHD")
    filtered <- .tapscan_exclude_similar(hits)
    expect_equal(filtered$Domain, c("G2-like_Domain", "WD40", "PHD"))
})


test_that(".tapscan_tokens() adds pseudo-domains", {
    hits <- data.frame(
        Gene = c("g1", "g1", "g2", "g3"),
        Domain = c("Myb_DNA-binding", "LIM", "Myb_DNA-binding", "LIM"),
        c_evalue = 1e-10, count = c(3, 2, 5, 1)
    )
    tokens <- .tapscan_tokens(hits)
    expect_equal(
        tokens$Token,
        c("MYB-3R", "Myb_DNA-binding", "two_or_more_LIM", "LIM",
          "Myb_DNA-binding", "LIM")
    )
    expect_equal(tokens$Gene, c("g1", "g1", "g1", "g1", "g2", "g3"))
})


test_that(".tapscan_tiebreak() chooses a single family", {
    expect_equal(.tapscan_tiebreak("A", ";d1;", "d1"), "A")

    # First token found only in the second candidate: switch
    expect_equal(
        .tapscan_tiebreak(c("A", "B"), c(";d1;", ";d2;d1;"), c("d2", "d1")),
        "B"
    )
    # First token found only in the first candidate: keep
    expect_equal(
        .tapscan_tiebreak(c("A", "B"), c(";d1;", ";d2;"), c("d1", "d2")),
        "A"
    )
    # Tokens shared by both candidates: keep
    expect_equal(
        .tapscan_tiebreak(c("A", "B"), c(";d1;", ";d1;"), "d1"), "A"
    )
})


test_that(".tapscan_family_map() merges families and adds subfamilies", {
    fams <- .tapscan_family_map(c("bZIP2", "MYB-2R", "NAC", "HRT"))
    expect_equal(fams$Family, c("bZIP", "MYB", "NAC", "ET"))
    expect_equal(fams$Subfamily, c("bZIP", "MYB-2R", "NAC", "ET"))
})


test_that(".classify_tapscan() handles input without classified genes", {
    annot <- data.frame(
        Gene = "g1", Domain = "AP2", qlen = 100, c_evalue = 1e-10,
        hmm_from = 1, hmm_to = 2
    )
    fams <- .classify_tapscan(annot)
    expect_equal(nrow(fams), 0)
    expect_equal(names(fams), c("Gene", "Family", "Subfamily"))
})


test_that(".add_tap_class() adds TAP classes", {
    fams <- data.frame(
        TAPscan_family = c("MYB", "HAT", "NAC", "PcG_FIE"),
        TAPscan_subfamily = c("MYB-2R", "GNAT", "NAC", "PcG_FIE")
    )
    res <- .add_tap_class(fams)
    expect_equal(res$TAP_class, c("TF", "TR", "TF", NA))
})


test_that("tap_class data set has classes for TAPscan families", {
    data(tap_class)
    expect_equal(
        names(tap_class),
        c("TAP_family", "TAP_class", "Description", "References")
    )
    expect_true(all(tap_class$TAP_class %in% c("TF", "TR", "PT")))
    expect_false(any(duplicated(tap_class$TAP_family)))
    expect_equal(tap_class[, 1:2], tap_class_map)

    # Most families produced by TAPscan rules have a TAP class
    fams <- unique(.tapscan_family_map(unique(tapscan_rules$family)))
    names(fams) <- c("TAPscan_family", "TAPscan_subfamily")
    fams <- .add_tap_class(fams)
    expect_true(mean(!is.na(fams$TAP_class)) > 0.9)
})
