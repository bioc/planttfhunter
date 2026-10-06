#----Load data------------------------------------------------------------------
# Helper function to create input for .classify_planttfclass()
make_families <- function(ptfdb, tap_fam, tap_sub = tap_fam) {
    data.frame(
        Gene = sprintf("g%d", seq_along(tap_fam)),
        PlantTFDB_subfamily = ptfdb,
        TAPscan_family = tap_fam,
        TAPscan_subfamily = tap_sub
    )
}

#----Start tests----------------------------------------------------------------
test_that("correspondence table has one superclass and class per family", {
    map <- unique(planttfclass_map[, c("family", "class", "superclass")])
    expect_false(any(duplicated(map$family)))
    expect_false(any(is.na(planttfclass_map$family)))
    expect_false(any(is.na(planttfclass_map$class)))
    expect_false(any(is.na(planttfclass_map$superclass)))

    data(planttfclass_scheme)
    expect_equal(
        planttfclass_scheme,
        setNames(planttfclass_map[, c(1, 2, 5, 4, 3, 6)],
                 names(planttfclass_scheme))
    )
})


test_that(".classify_planttfclass() adds Plant-TFClass classifications", {
    fams <- make_families(
        ptfdb = c("WRKY", NA, NA, NA, NA, NA),
        tap_fam = c("WRKY", "bHLH", "bHLH", "HDZ", "LFY", "FHA"),
        tap_sub = c("WRKY", "bHLH", "bHLH_TCP", "C3HDZ", "LFY", "FHA")
    )
    res <- .classify_planttfclass(fams)

    expect_equal(
        names(res),
        c(names(fams), "PlantTFClass_superclass", "PlantTFClass_class",
          "PlantTFClass_family")
    )
    # Classes without families (e.g., LEAFY) have the class name as family
    expect_equal(
        res$PlantTFClass_family,
        c("WRKY", "bHLH", "TCP", "HD-ZIP", "LEAFY", NA)
    )
    expect_equal(
        res$PlantTFClass_class,
        c("GCM domain factors", "Basic helix-loop-helix factors (bHLH)",
          "TCP", "Homeo domain factors", "LEAFY", NA)
    )
    expect_equal(res$PlantTFClass_superclass[3], "Beta-sheet binding to DNA")
    expect_true(is.na(res$PlantTFClass_superclass[6]))
})


test_that("subfamily-specific matches take precedence over family matches", {
    fams <- make_families(
        ptfdb = NA,
        tap_fam = c("MYB", "MYB", "C2H2", "C2H2"),
        tap_sub = c("MYB-2R", "MYB-related", "C2H2", "C2H2_IDD")
    )
    res <- .classify_planttfclass(fams)
    expect_equal(
        res$PlantTFClass_family, c("MYB", "MYB-related", "C2H2", "IDD")
    )
})


test_that("ambiguous classifications are resolved with PlantTFDB families", {
    fams <- make_families(
        ptfdb = c("ERF", "AP2", "RAV", "M-type", "MIKC", "B3", "RAV",
                  "G2-like;ERF"),
        tap_fam = c("AP2", "AP2", "AP2", "MADS", "MADS", "ABI3/VP1",
                    "ABI3/VP1", "AP2"),
        tap_sub = c(
            "AP2", "AP2", "CRF", "MADS", "MADS", "ABI3/VP1", "ABI3/VP1", "AP2"
        )
    )
    res <- .classify_planttfclass(fams)
    expect_equal(
        res$PlantTFClass_family,
        c("ERF/DREB", "AP2", "RAV", "Type I", "Type II", "LAV", "RAV",
          "ERF/DREB")
    )
    expect_equal(res$PlantTFClass_class[3], "B3")
    expect_equal(
        res$PlantTFClass_superclass[3], "Beta-barrel DNA-binding domains"
    )
})


test_that("unresolved ambiguous classifications are separated by ';'", {
    fams <- make_families(
        ptfdb = c(NA, "bZIP", NA), tap_fam = c("AP2", "AP2", "MADS"),
        tap_sub = c("AP2", "AP2", "MADS")
    )
    res <- .classify_planttfclass(fams)

    expect_equal(res$PlantTFClass_family[1], "AP2;ERF/DREB;RAV")
    expect_equal(res$PlantTFClass_family[2], "AP2;ERF/DREB;RAV")
    expect_equal(res$PlantTFClass_class[1], "AP2/EREBP;B3")
    expect_equal(
        res$PlantTFClass_superclass[1],
        paste0(
            "Beta-hairpin exposed by an alpha/beta-scaffold;",
            "Beta-barrel DNA-binding domains"
        )
    )
    expect_equal(res$PlantTFClass_family[3], "Type I;Type II")
    expect_equal(res$PlantTFClass_class[3], "MADS box factors")
})


test_that(".resolve_planttfclass() handles empty input", {
    res <- .resolve_planttfclass(planttfclass_map[0, ], NA)
    expect_equal(res, rep(NA_character_, 3))

    fams <- make_families(character(0), character(0), character(0))
    expect_equal(nrow(.classify_planttfclass(fams)), 0)
})
