
#----Load data------------------------------------------------------------------
data(gsu_annotation)
domain_annotation <- gsu_annotation

#----Start tests----------------------------------------------------------------
test_that("classify_tfs() returns a 2-column data frame", {
  families <- classify_tfs(domain_annotation)
  expect_equal(class(families), "data.frame")
  expect_equal(ncol(families), 2)
})

test_that("list_domains() returns a named character vector", {
    dom_dbd <- list_domains(show = "dbd")
    dom_forbidden <- list_domains(show = "forbidden")
    dom_aux <- list_domains(show = "auxiliary")
    
    expect_equal(length(dom_dbd), 47)
    expect_equal(length(dom_forbidden), 5)
    expect_equal(length(dom_aux), 10)
    expect_equal(class(dom_dbd), "character")
    expect_equal(class(dom_forbidden), "character")
    expect_equal(class(dom_aux), "character")
})


test_that("check_ap2_erf() returns domain classification", {
    ct <- check_ap2_erf("PF00847")
    cf <- check_ap2_erf("fakedomain")
    
    expect_equal(ct, "ERF")
    expect_equal(cf, NULL)
})

test_that("check_b3() returns domain classification", {
    ct <- check_b3("PF02362")
    cf <- check_b3("fakedomain")
    
    expect_equal(ct, "B3")
    expect_equal(cf, NULL)
})


test_that("check_c2c2() returns domain classification", {
    ct <- check_c2c2("PF00320")
    cf <- check_c2c2("fakedomain")
    
    expect_equal(ct, "GATA")
    expect_equal(cf, NULL)
})


test_that("check_garp() returns domain classification", {
    ct <- check_garp("G2-like")
    cf <- check_garp("fakedomain")
    
    expect_equal(ct, "G2-like")
    expect_equal(cf, NULL)
})


test_that("check_hb() returns domain classification", {
    ct <- check_hb("PF00046")
    cf <- check_hb("fakedomain")
    
    expect_equal(ct, "HB-other")
    expect_equal(cf, NULL)
})


test_that("check_mads() returns domain classification", {
    ct <- check_mads("PF00319")
    cf <- check_mads("fakedomain")
    
    expect_equal(ct, "M-type")
    expect_equal(cf, NULL)
})


test_that("check_myb() returns domain classification", {
    ct <- check_myb("PF00249")
    cf <- check_myb("fakedomain")
    
    expect_equal(ct, "MYB-related")
    expect_equal(cf, NULL)
})


test_that("check_nf_y() returns domain classification", {
    ct <- check_nf_y("PF02045")
    cf <- check_nf_y("fakedomain")
    
    expect_equal(ct, "NF-YA")
    expect_equal(cf, NULL)
})


test_that("check_smallfams() returns domain classification", {
    ct <- check_smallfams("PF00010")
    cf <- check_smallfams("fakedomain")
    
    expect_equal(ct, "bHLH")
    expect_equal(cf, NULL)
})
