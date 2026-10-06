#----Load data------------------------------------------------------------------
data(gsu)
data(gsu_annotation)
seq <- gsu[1:5]

#----Start tests----------------------------------------------------------------
test_that("annotate_domains() returns a list of domain annotations", {
    skip_if_not(hmmer_is_installed())

    annot <- annotate_domains(seq, threads = 2)
    cols <- c("Gene", "Domain", "qlen", "c_evalue", "hmm_from", "hmm_to")
    expect_equal(names(annot), c("PlantTFDB", "TAPscan"))
    expect_equal(names(annot$PlantTFDB), cols)
    expect_equal(names(annot$TAPscan), cols)
    expect_true(nrow(annot$PlantTFDB) >= 2)
    expect_true(nrow(annot$TAPscan) > 0)

    ref <- gsu_annotation$TAPscan[gsu_annotation$TAPscan$Gene %in% names(seq), ]
    key <- function(x) sort(paste(x$Gene, x$Domain, x$hmm_from, x$hmm_to))
    expect_equal(key(annot$TAPscan), key(ref))

    expect_error(annotate_domains(as.character(seq)), "AAStringSet")
})


test_that("annotate_domains() fails if HMMER is not installed", {
    path <- Sys.getenv("PATH")
    Sys.setenv(PATH = "")
    installed <- hmmer_is_installed()
    err <- tryCatch(annotate_domains(seq), error = conditionMessage)
    Sys.setenv(PATH = path)

    expect_false(installed)
    expect_match(err, "Could not find HMMER")
})
