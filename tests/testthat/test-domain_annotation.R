
#----Load data------------------------------------------------------------------
data(gsu)
seq <- gsu[1:2]

#----Start tests----------------------------------------------------------------
test_that("annotate_pfam() returns a data frame of domain annotation", {
    annot <- annotate_pfam(seq)

    annot_local <- data.frame(
        Gene = c("gsu04730.1", "gsu00340.1"),
        Domain = c("PF02042", "PF00076")
    )
    if(hmmer_is_installed()) {
        annot_local <- annotate_pfam(seq, mode = "local")
    }
    
    expect_equal(class(annot), "data.frame")
    expect_equal(ncol(annot), 2)
    expect_equal(names(annot), c("Gene", "Domain"))
    expect_true(nrow(annot) >= 2)
    
    expect_equal(class(annot_local), "data.frame")
    expect_equal(ncol(annot_local), 2)
    expect_equal(names(annot_local), c("Gene", "Domain"))
    expect_true(nrow(annot_local) >= 2)
    
    expect_error(annotate_pfam(seq, mode = "native"))
})
