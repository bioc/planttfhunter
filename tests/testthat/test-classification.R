
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
    expect_error(list_domains(show = "anything"))
})


test_that("check_ap2_erf() returns domain classification", {
    ct <- check_ap2_erf("PF00847")
    ct2 <- check_ap2_erf(rep("PF00847", 2))
    ct3 <- check_ap2_erf(c("PF02362", "PF00847"))
    cf <- check_ap2_erf("fakedomain")
    
    expect_equal(ct, "ERF")
    expect_equal(ct2, "AP2")
    expect_equal(ct3, "RAV")
    expect_equal(cf, NULL)
})


test_that("check_b3() returns domain classification", {
    ct <- check_b3("PF02362")
    ct2 <- check_b3(c("PF06507", "PF02362"))
    cf <- check_b3("fakedomain")
    
    expect_equal(ct, "B3")
    expect_equal(ct2, "ARF")
    expect_equal(cf, NULL)
})


test_that("check_c2c2() returns domain classification", {
    ct <- check_c2c2("PF00320")
    ct2 <- check_c2c2(c("PF00643", "PF06203"))
    ct3 <- check_c2c2("PF02701")
    ct4 <- check_c2c2("PF06943")
    ct5 <- check_c2c2("PF04690")
    cf <- check_c2c2("fakedomain")
    
    expect_equal(ct, "GATA")
    expect_equal(ct2, "CO-like")
    expect_equal(ct3, "Dof")
    expect_equal(ct4, "LSD")
    expect_equal(ct5, "YABBY")
    expect_equal(cf, NULL)
})


test_that("check_garp() returns domain classification", {
    ct <- check_garp("G2-like")
    ct2 <- check_garp(c("PF00072", "G2-like"))
    cf <- check_garp("fakedomain")
    
    expect_equal(ct, "G2-like")
    expect_equal(ct2, "ARR-B")
    expect_equal(cf, NULL)
})


test_that("check_hb() returns domain classification", {
    ct <- check_hb("PF00046")
    ct2 <- check_hb(c("PF00046", "HD-ZIP_I/II"))
    ct3 <- check_hb(c("PF03789", "PF00046"))
    ct4 <- check_hb(c("PF00628", "PF00046"))
    ct5 <- check_hb(c("BELL", "PF00046"))
    ct6 <- check_hb(c("Wus_type_Homeobox", "PF00046"))
    ct7 <- check_hb(c("PF01852", "PF00046"))
    cf <- check_hb("fakedomain")
    
    expect_equal(ct, "HB-other")
    expect_equal(ct2, "HD-ZIP")
    expect_equal(ct3, "TALE")
    expect_equal(ct4, "HB-PHD")
    expect_equal(ct5, "TALE")
    expect_equal(ct6, "WOX")
    expect_equal(ct7, "HD-ZIP")
    expect_equal(cf, NULL)
})


test_that("check_mads() returns domain classification", {
    ct <- check_mads("PF00319")
    ct2 <- check_mads(c("PF01486", "PF00319"))
    cf <- check_mads("fakedomain")
    
    expect_equal(ct, "M-type")
    expect_equal(ct2, "MIKC")
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
    ct2 <- check_smallfams("PF06217")
    ct3 <- check_smallfams("PF05687")
    ct4 <- check_smallfams(c("PF00642", "PF00271"))
    ct5 <- check_smallfams("PF03859")
    ct6 <- check_smallfams(rep("PF00643", 2))
    ct7 <- check_smallfams("PF04873")
    ct8 <- check_smallfams("PF03101")
    ct9 <- check_smallfams("PF04504")
    ct10 <- check_smallfams("PF03514")
    ct11 <- check_smallfams(c("PF08879", "PF08880"))
    ct12 <- check_smallfams("HRT-like")
    ct13 <- check_smallfams("PF03195")
    ct14 <- check_smallfams("PF01698")
    ct15 <- check_smallfams("PF02365")
    ct16 <- check_smallfams("PF08744")
    ct17 <- check_smallfams("PF04689")
    ct18 <- check_smallfams("SAP")
    ct19 <- check_smallfams("PF03110")
    ct20 <- check_smallfams("PF05142")
    ct21 <- check_smallfams("STAT")
    ct22 <- check_smallfams("PF03634")
    ct23 <- check_smallfams("trihelix")
    ct24 <- check_smallfams("VOZ")
    ct25 <- check_smallfams("PF08536")
    ct26 <- check_smallfams("PF03106")
    ct27 <- check_smallfams("PF04770")
    cf <- check_smallfams("fakedomain")
    
    expect_equal(ct, "bHLH")
    expect_equal(ct2, "BBR-BPC")
    expect_equal(ct3, "BES1")
    expect_equal(ct4, "C3H")
    expect_equal(ct5, "CAMTA")
    expect_equal(ct6, "DBB")
    expect_equal(ct7, "EIL")
    expect_equal(ct8, "FAR1")
    expect_equal(ct9, "GeBP")
    expect_equal(ct10, "GRAS")
    expect_equal(ct11, "GRF")
    expect_equal(ct12, "HRT-like")
    expect_equal(ct13, "LBD")
    expect_equal(ct14, "LFY")
    expect_equal(ct15, "NAC")
    expect_equal(ct16, "NZZ/SPL")
    expect_equal(ct17, "S1Fa-like")
    expect_equal(ct18, "SAP")
    expect_equal(ct19, "SBP")
    expect_equal(ct20, "SRS")
    expect_equal(ct21, "STAT")
    expect_equal(ct22, "TCP")
    expect_equal(ct23, "Trihelix")
    expect_equal(ct24, "VOZ")
    expect_equal(ct25, "Whirly")
    expect_equal(ct26, "WRKY")
    expect_equal(ct27, "ZF-HD")

    expect_equal(cf, NULL)
})
