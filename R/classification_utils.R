

#' List possible domains used for TF classification
#'
#' @param show Character specifying what to show. One of 'dbd', 'auxiliary',
#' 'forbidden', or 'all'. Default: 'dbd'.
#'
#' @return A named character vector of domains and their IDs.
#' @noRd
list_domains <- function(show = "dbd") {
  domains <- c(
    AP2 = "PF00847", B3 = "PF02362", auxin_resp = "PF06507", 
    GAGA_bind = "PF06217", DUF822 = "PF05687", HLH = "PF00010",
    bZIP_1 = "PF00170", zf_B_box = "PF00643", CCT = "PF06203", 
    Zf_Dof = "PF02701", GATA_zf = "PF00320", Zf_LSD1 = "PF06943",
    Peptidase_C14 = "PF00656", YABBY = "PF04690", Zf_C2H2 =  "PF00096",
    RNase_T = "PF00929", Zf_CCCH =  "PF00642", RRM_1 = "PF00076",
    Helicase_C = "PF00271", CG1 = "PF03859", TCR = "PF03638",
    E2F_TDP = "PF02319", EIN3 = "PF04873", FAR1 = "PF03101",
    STAT = "PF02864", Response_reg = "PF00072", DUF573 = "PF04504",
    GRAS = "PF03514", WRC = "PF08879", QLQ = "PF08880",
    Homeobox = "PF00046", SMART = "PF01852", ELK = "PF03789",
    PHD = "PF00628", HSF = "PF00447", DUF260 = "PF03195",
    FLO_LFY = "PF01698", SRF_TF = "PF00319", K_box = "PF01486",
    MYB_DNAbind = "PF00249", SWIRM = "PF04433", NAM = "PF02365",
    Zf_NF_X1 = "PF01422", CBFB_NFYA = "PF02045", RWP_RK = "PF02042",
    NOZZLE = "PF08744", S1FA = "PF04689", SBP = "PF03110",
    DUF702 = "PF05142", TCP = "PF03634", WRKY = "PF03106",
    Whirly = "PF08536", ZF_HD_dimer = "PF04770",
    G2_like = "G2-like", HD_ZIP = "HD-ZIP",
    HRT_like = "HRT-like", NF_YB = "NF-YB", NF_YC = "NF-YC",
    SAP = "SAP", Trihelix = "Trihelix", VOZ = "VOZ", WOX = "WOX"
  )
  aux_domains <- c("PF06507", "PF06203", "PF00072", "PF08880", "PF01852", 
                   "HD-ZIP", "PF03789", "WOX", "PF00628", "PF01486")
  forb_domains <- c("PF04433", "PF00271", "PF00076", "PF00929", "PF00656")
  
  if(show == "forbidden") {
    fdomains <- domains[domains %in% forb_domains]
  } else if(show == "auxiliary") {
    fdomains <- domains[domains %in% aux_domains]
  } else if(show == "dbd") {
    fdomains <- domains[!(domains %in% c(forb_domains, aux_domains))]
  } else {
    stop("Invalid argument to parameter show.")
  }
  return(fdomains)
}


#' Check if TF belongs to the AP2/ERF family
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_ap2_erf <- function(domains = NULL) {
  len_ap2 <- length(domains[domains == "PF00847"])
  len_b3 <- length(domains[domains == "PF02362"])
  
  if(len_ap2 >= 2) {
    fam <- "AP2"
  } else if(len_ap2 == 1 & len_b3 == 0) {
    fam <- "ERF"
  } else if(len_ap2 != 0 & len_b3 != 0) {
    fam <- "RAV"
  } else {
    fam <- NULL
  }
  return(fam)
}


#' Check if TF belongs to the B3 family
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_b3 <- function(domains = NULL) {
  len_ap2 <- length(domains[domains == "PF00847"])
  len_auxin <- length(domains[domains == "PF06507"])
  len_b3 <- length(domains[domains == "PF02362"])
  
  if(len_b3 != 0 & len_auxin != 0) {
    fam <- "ARF"
  } else if(len_b3 != 0 & len_ap2 == 0 & len_auxin == 0) {
    fam <- "B3"
  } else {
    fam <- NULL
  }
  return(fam)
}


#' Check if TF belongs to the C2C2 family
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_c2c2 <- function(domains = NULL) {
  len_zf_bbox <- length(domains[domains == "PF00643"])
  len_cct <- length(domains[domains == "PF06203"])
  len_zf_dof <- length(domains[domains == "PF02701"])
  len_gata_zf <- length(domains[domains == "PF00320"])
  len_zf_lsd1 <- length(domains[domains == "PF06943"])
  len_pep_c14 <- length(domains[domains == "PF00656"])
  len_yabby <- length(domains[domains == "PF04690"])
  
  if(len_zf_bbox != 0 & len_cct != 0) {
    fam <- "CO-like"
  } else if(len_zf_dof != 0) {
    fam <- "Dof"
  } else if(len_gata_zf != 0) {
    fam <- "GATA"
  } else if(len_zf_lsd1 != 0 & len_pep_c14 == 0) {
    fam <- "LSD"
  } else if(len_yabby != 0) {
    fam <- "YABBY"
  } else {
    fam <- NULL
  }
  
  return(fam)
}


#' Check if TF belongs to the GARP family
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_garp <- function(domains = NULL) {
    len_g2like <- length(domains[domains == "G2-like"])
    len_resp_rep <- length(domains[domains == "PF00072"])
    
    if(len_g2like != 0 & len_resp_rep != 0) {
        fam <- "ARR-B"
    } else if(len_g2like != 0 & len_resp_rep == 0) {
        fam <- "G2-like"
    } else {
        fam <- NULL
    }
    
    return(fam)
}


#' Check if TF belongs to the HB family
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_hb <- function(domains = NULL) {
    len_homeobox <- length(domains[domains == "PF00046"])
    len_hdzip <- length(domains[domains == "HD-ZIP_I/II"])
    len_smart <- length(domains[domains == "PF01852"])
    len_elk <- length(domains[domains == "PF03789"])
    len_bell <- length(domains[domains == "BELL"])
    len_wox <- length(domains[domains == "Wus_type_Homeobox"])
    len_phd <- length(domains[domains == "PF00628"])
    only_homeobox <- unique(domains)
    only_homeobox <- length(only_homeobox) == 1 & "PF00046" %in% only_homeobox
    
    if(len_homeobox != 0 & len_hdzip != 0) {
        fam <- "HD-ZIP"
    } else if(len_homeobox != 0 & len_smart != 0) {
        fam <- "HD-ZIP"
    } else if(len_homeobox != 0 & len_elk != 0) {
        fam <- "TALE"
    } else if(len_homeobox != 0 & len_bell != 0) {
        fam <- "TALE"
    } else if(len_homeobox != 0 & len_wox != 0) {
        fam <- "WOX"
    } else if(len_homeobox != 0 & len_phd != 0) {
        fam <- "HB-PHD"
    } else if(only_homeobox) {
        fam <- "HB-other"
    } else {
        fam <- NULL
    }
    
    return(fam)
}


#' Check if TF belongs to the MADS-box family
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_mads <- function(domains = NULL) {
    len_srf_tf <- length(domains[domains == "PF00319"])
    len_kbox <- length(domains[domains == "PF01486"])
    
    if(len_srf_tf != 0 & len_kbox == 0) {
        fam <- "M-type"
    } else if(len_srf_tf != 0 & len_kbox != 0) {
        fam <- "MIKC"
    } else {
        fam <- NULL
    }
    
    return(fam)
}


#' Check if TF belongs to the MYB family
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_myb <- function(domains = NULL) {
    len_myb <- length(domains[domains == "PF00249"])
    len_swirm <- length(domains[domains == "PF04433"])
    
    if(len_myb == 1 & len_swirm == 0) {
        fam <- "MYB-related"
    } else if(len_myb > 1 & len_swirm == 0) {
        fam <- "MYB"
    } else {
        fam <- NULL
    }
    
    return(fam)
}


#' Check if TF belongs to the NF-Y family
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_nf_y <- function(domains = NULL) {
    len_cbfb_nfya <- length(domains[domains == "PF02045"])
    len_nf_yb <- length(domains[domains == "NF-YB"])
    len_nf_yc <- length(domains[domains == "NF-YC"])
    
    if(len_cbfb_nfya != 0) {
        fam <- "NF-YA"
    } else if(len_nf_yb != 0) {
        fam <- "NF-YB"
    } else if(len_nf_yc != 0) {
        fam <- "NF-YC"
    } else {
        fam <- NULL
    }
    
    return(fam)
}


#' Check if TF belongs to one of the smaller families
#'
#' @param domains Character scalar or vector of domains for a given gene.
#' 
#' @return A character scalar with the TF family or NULL if it belongs to a
#' different family.
#' @noRd
check_smallfams <- function(domains = NULL) {
    count <- function(x) { return(length(domains[domains == x])) }
    
    if(count("PF06217") != 0) { # GAGA_bind
        fam <- "BBR-BPC"
    } else if(count("PF05687") != 0) { # DUF822
        fam <- "BES1"
    } else if(count("PF00010") != 0) { # HLH
        fam <- "bHLH"
    } else if(count("PF00170") != 0) { # bZIP_1
        fam <- "bZIP"
    } else if(count("PF00096") != 0 & count("PF00929") == 0) { # zf-C2H2, RNase_T
        fam <- "C2H2"
    } else if(count("PF00642") != 0 & count("PF00076") == 0) { # zf-CCCH, RRM_1
        fam <- "C3H"
    } else if(count("PF00642") != 0 & count("PF00271") == 0) { # ., Helicase_C
        fam <- "C3H"
    } else if(count("PF03859") != 0) { # CG1
        fam <- "CAMTA"
    } else if(count("PF03638") != 0) { # TCR
        fam <- "CPP"
    } else if(count("PF00643") > 1) { # zf-B_box
        fam <- "DBB"
    } else if(count("PF02319") != 0) { # E2F_TDP
        fam <- "E2F/DP"
    } else if(count("PF04873") != 0) { # EIN3
        fam <- "EIL"
    } else if(count("PF03101") != 0) { # FAR1
        fam <- "FAR1"
    } else if(count("PF04504") != 0) { # DUF573
        fam <- "GeBP"
    } else if(count("PF03514") != 0) { # GRAS
        fam <- "GRAS"
    } else if(count("PF08879") != 0 & count("PF08880") != 0) { # WRC, QLQ
        fam <- "GRF"
    } else if(count("HRT-like") != 0) { # HRT-like, self-built
        fam <- "HRT-like"
    } else if(count("PF00447") != 0) { # HSF_dna_bind
        fam <- "HSF"
    } else if(count("PF03195") != 0) { # DUF260
        fam <- "LBD"
    } else if(count("PF01698") != 0) { # FLO_LFY
        fam <- "LFY"
    } else if(count("PF02365") != 0) { # NAM
        fam <- "NAC"
    } else if(count("PF01422") != 0) { # zf-NF-X1
        fam <- "NF-X1"
    } else if(count("PF02042") != 0) { # RWP-RK
        fam <- "Nin-like"
    } else if(count("PF08744") != 0) { # NOZZLE
        fam <- "NZZ/SPL"
    } else if(count("PF04689") != 0) { # S1FA
        fam <- "S1Fa-like"
    } else if(count("SAP") != 0) { # SAP
        fam <- "SAP"
    } else if(count("PF03110") != 0) { # SBP
        fam <- "SBP"
    } else if(count("PF05142") != 0) { # DUF702
        fam <- "SRS"
    } else if(count("STAT") != 0) { # STAT
        fam <- "STAT"
    } else if(count("PF03634") != 0) { # TCP
        fam <- "TCP"
    } else if(count("trihelix") != 0) { # Trihelix
        fam <- "Trihelix"
    } else if(count("VOZ") != 0) { # VOZ
        fam <- "VOZ"
    } else if(count("PF08536") != 0) { # Whirly
        fam <- "Whirly"
    } else if(count("PF03106") != 0) { # WRKY
        fam <- "WRKY"
    } else if(count("PF04770") != 0) { # ZF-HD_dimer
        fam <- "ZF-HD"
    } else {
        fam <- NULL
    }
    return(fam)
}



