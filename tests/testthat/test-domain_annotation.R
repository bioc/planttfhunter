
#----Load data------------------------------------------------------------------
data(gsu)
seq <- gsu[1:2]

#----Start tests----------------------------------------------------------------
test_that("annotate_pfam() returns a data frame of domain annotation", {
  annot <- annotate_pfam(seq)
  expect_equal(class(annot), "data.frame")
  expect_equal(ncol(annot), 2)
})
