
#' List domains used for TF classification with PlantTFDB's scheme
#'
#' @return A named character vector of domain IDs (or domain names, for
#' self-built domains), including DNA-binding, auxiliary, and forbidden
#' domains (see \code{planttfdb_scheme}).
#' @noRd
.list_domains <- function() {

    domains <- c(
        AP2 = "PF00847", B3 = "PF02362", auxin_resp = "PF06507",
        GAGA_bind = "PF06217", DUF822 = "PF05687", HLH = "PF00010",
        bZIP_1 = "PF00170", zf_B_box = "PF00643", CCT = "PF06203",
        Zf_Dof = "PF02701", GATA_zf = "PF00320", Zf_LSD1 = "PF06943",
        Peptidase_C14 = "PF00656", YABBY = "PF04690", Zf_C2H2 = "PF00096",
        RNase_T = "PF00929", Zf_CCCH = "PF00642", RRM_1 = "PF00076",
        Helicase_C = "PF00271", CG1 = "PF03859", TCR = "PF03638",
        E2F_TDP = "PF02319", EIN3 = "PF04873", FAR1 = "PF03101",
        Response_reg = "PF00072", DUF573 = "PF04504",
        GRAS = "PF03514", WRC = "PF08879", QLQ = "PF08880",
        Homeobox = "PF00046", SMART = "PF01852", ELK = "PF03789",
        PHD = "PF00628", HSF = "PF00447", DUF260 = "PF03195",
        FLO_LFY = "PF01698", SRF_TF = "PF00319", K_box = "PF01486",
        MYB_DNAbind = "PF00249", SWIRM = "PF04433", NAM = "PF02365",
        Zf_NF_X1 = "PF01422", CBFB_NFYA = "PF02045", RWP_RK = "PF02042",
        NOZZLE = "PF08744", S1FA = "PF04689", SBP = "PF03110",
        DUF702 = "PF05142", TCP = "PF03634", WRKY = "PF03106",
        Whirly = "PF08536", ZF_HD_dimer = "PF04770",
        G2_like = "G2-like", HD_ZIP = "HD-ZIP_I/II", BELL = "BELL",
        WOX = "Wus_type_Homeobox", HRT_like = "HRT-like",
        NF_YB = "NF-YB", NF_YC = "NF-YC", SAP = "SAP", STAT = "STAT",
        Trihelix = "trihelix", VOZ = "VOZ"
    )
    
    return(domains)
}


#' Get rules for TF classification using PlantTFDB's scheme
#'
#' @param m A gene x domain matrix of domain counts.
#'
#' @return A list of named lists. Each element of the outer list represents
#' a group of mutually exclusive families, and each inner list contains
#' logical vectors (one per family, in order of priority) indicating whether
#' each gene (row of \strong{m}) satisfies the rule for that family.
#' Within a group, genes are assigned to the first family whose rule they
#' satisfy.
#' @noRd
.planttfdb_rules <- function(m) {

    n <- function(d) m[, d]
    has <- function(d) m[, d] > 0
    hb <- has("PF00046")
    only_hb <- hb & rowSums(m > 0) == 1

    rules <- list(
        ap2_erf = list(
            AP2 = n("PF00847") >= 2,
            ERF = n("PF00847") == 1 & !has("PF02362"),
            RAV = has("PF00847") & has("PF02362")
        ),
        b3 = list(
            ARF = has("PF02362") & has("PF06507"),
            B3 = has("PF02362") & !has("PF00847") & !has("PF06507")
        ),
        c2c2 = list(
            `CO-like` = has("PF00643") & has("PF06203"),
            Dof = has("PF02701"),
            GATA = has("PF00320"),
            LSD = has("PF06943") & !has("PF00656"),
            YABBY = has("PF04690")
        ),
        garp = list(
            `ARR-B` = has("G2-like") & has("PF00072"),
            `G2-like` = has("G2-like") & !has("PF00072")
        ),
        hb = list(
            `HD-ZIP` = hb & (has("HD-ZIP_I/II") | has("PF01852")),
            TALE = hb & (has("PF03789") | has("BELL")),
            WOX = hb & has("Wus_type_Homeobox"),
            `HB-PHD` = hb & has("PF00628"),
            `HB-other` = only_hb
        ),
        mads = list(
            `M-type` = has("PF00319") & !has("PF01486"),
            MIKC = has("PF00319") & has("PF01486")
        ),
        myb = list(
            `MYB-related` = n("PF00249") == 1 & !has("PF04433"),
            MYB = n("PF00249") > 1 & !has("PF04433")
        ),
        nf_y = list(
            `NF-YA` = has("PF02045"),
            `NF-YB` = has("NF-YB"),
            `NF-YC` = has("NF-YC")
        )
    )
    rules$smallfams <- .planttfdb_smallfam_rules(n, has)
    
    return(rules)
}


#' Get rules for smaller TF families using PlantTFDB's scheme
#'
#' @param n Function that takes a domain ID and returns the number of
#' times the domain was found in each gene.
#' @param has Function that takes a domain ID and returns a logical vector
#' indicating whether each gene has the domain.
#'
#' @return A named list of logical vectors (one per family, in order of
#' priority) as in \code{.planttfdb_rules()}.
#' @noRd
.planttfdb_smallfam_rules <- function(n, has) {

    rules <- list(
        `BBR-BPC` = has("PF06217"),
        BES1 = has("PF05687"),
        bHLH = has("PF00010"),
        bZIP = has("PF00170"),
        C2H2 = has("PF00096") & !has("PF00929"),
        C3H = has("PF00642") & !has("PF00076") & !has("PF00271"),
        CAMTA = has("PF03859"),
        CPP = has("PF03638"),
        DBB = n("PF00643") > 1 & !has("PF06203"),
        `E2F/DP` = has("PF02319"),
        EIL = has("PF04873"),
        FAR1 = has("PF03101"),
        GeBP = has("PF04504"),
        GRAS = has("PF03514"),
        GRF = has("PF08879") & has("PF08880"),
        `HRT-like` = has("HRT-like"),
        HSF = has("PF00447"),
        LBD = has("PF03195"),
        LFY = has("PF01698"),
        NAC = has("PF02365"),
        `NF-X1` = has("PF01422"),
        `Nin-like` = has("PF02042"),
        `NZZ/SPL` = has("PF08744"),
        `S1Fa-like` = has("PF04689"),
        SAP = has("SAP"),
        SBP = has("PF03110"),
        SRS = has("PF05142"),
        STAT = has("STAT"),
        TCP = has("PF03634"),
        Trihelix = has("trihelix"),
        VOZ = has("VOZ"),
        Whirly = has("PF08536"),
        WRKY = has("PF03106"),
        `ZF-HD` = has("PF04770")
    )
    
    return(rules)
}

#' Assign each gene to the first family whose rule it satisfies
#'
#' @param rules A named list of logical vectors (one per family, in order
#' of priority).
#'
#' @return A character vector with the family of each gene, or NA if a gene
#' does not satisfy any rule.
#' @noRd
.first_match <- function(rules) {

    init <- rep(NA_character_, length(rules[[1]]))
    fam <- Reduce(function(out, f) {
        out[is.na(out) & rules[[f]]] <- f
        return(out)
    }, names(rules), init)
    
    return(fam)
}


#' Identify TFs and classify them in families using PlantTFDB's scheme
#'
#' @param domain_annotation A data frame with variables
#' \strong{Gene} and \strong{Domain}, as in the element \strong{PlantTFDB}
#' of the output of \code{annotate_domains()}.
#'
#' @return A 3-column data frame with variables \strong{Gene},
#' \strong{Family}, and \strong{Subfamily}. As in PlantTFDB's classification
#' scheme, families without subfamilies have the family name as subfamily.
#' Genes that are assigned to more than one family have families (and
#' subfamilies, in the same order) separated by ';'.
#' @noRd
.classify_planttfdb <- function(domain_annotation = NULL) {

    dom <- .list_domains()
    annot <- domain_annotation[domain_annotation$Domain %in% dom, ]
    if(nrow(annot) == 0) {
        return(data.frame(
            Gene = character(0), Family = character(0),
            Subfamily = character(0)
        ))
    }

    # Gene x domain matrix of domain counts
    m <- unclass(table(
        factor(annot$Gene), factor(annot$Domain, levels = unname(dom))
    ))

    # Assign one subfamily per group, get their families, and then combine
    subfams <- lapply(.planttfdb_rules(m), .first_match)
    fams <- lapply(subfams, function(x) {
        planttfdb_families$Family[match(x, planttfdb_families$Subfamily)]
    })
    combine <- function(x) {
        Reduce(function(a, b) {
            ifelse(is.na(a), b, ifelse(is.na(b), a, paste(a, b, sep = ";")))
        }, x)
    }

    fam_df <- data.frame(
        Gene = rownames(m), Family = combine(fams),
        Subfamily = combine(subfams)
    )
    fam_df <- fam_df[!is.na(fam_df$Family), ]
    rownames(fam_df) <- NULL
    
    return(fam_df)
}
