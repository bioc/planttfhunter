
#----Load data------------------------------------------------------------------
data(gsu_annotation)
domain_annotation <- gsu_annotation

#----Start tests----------------------------------------------------------------
test_that("classify_tfs() returns a 2-column data frame", {
  families <- classify_tfs(domain_annotation)
  expect_equal(class(families), "data.frame")
  expect_equal(ncol(families), 2)
})
