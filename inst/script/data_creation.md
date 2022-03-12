Data creation
================

# Data in data/

## gsu.rda

This object contains the protein sequences of the algae *Galdieria
sulphuraria*, and it was downloaded from PLAZA Diatoms.

``` r
gsu <- Biostrings::readAAStringSet(
  "ftp://ftp.psb.ugent.be/pub/plaza/plaza_diatoms_01/Fasta/proteome.selected_transcript.gsu.fasta.gz"
)
names(gsu) <- gsub("\ .*", "", names(gsu))

usethis::use_data(
  gsu, compress = "xz"
)
```

## gsu_pfam.rda

``` r
data(gsu)
gsu_annotation <- annotate_pfam(gsu, mode = "local", evalue = 1e-05)

usethis::use_data(
  gsu_annotation, compress = "xz", overwrite = TRUE
)
```

## classification_scheme.rda

The classification scheme is the same as implemented in PlantTFDB.

``` r
library(tidyverse)
table <- rvest::read_html("http://planttfdb.gao-lab.org/help_famschema.php") %>%
  rvest::html_table()
table <- as.data.frame(table[[3]])

table <- table[-1, -1]
names(table) <- c("Family", "Subfamily", "DBD", "Auxiliary", "Forbidden")
table[table == ""] <- NA

table$DBD <- gsub("\\r\\n\\t\\t\\r\\n\\t\\t", " or ", table$DBD)
table$DBD <- gsub("\ \\(self-build\\)", "", table$DBD)
table$Auxiliary <- gsub("\ \\(self-build\\)", "", table$Auxiliary)
table$Auxiliary <- gsub("BELLor", "BELL or", table$Auxiliary)
table$Forbidden <- gsub("\\r\\n\\t\\tor", " or ", table$DBD)

classification_scheme <- table

usethis::use_data(classification_scheme,
                  compress = "xz")
```

# Data in inst/extdata

## \*.hmm files (profile HMMs)

The seed alignments for all TF families were downloaded from PFAM and
converted to profile HMMs using the command `hmmbuild` from HMMER.

``` r
library(tidyverse)
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
  Whirly = "PF08536", ZF_HD_dimer = "PF04770"
)

# Links to alignments of self-built domains in PlantTFDB (don't have PFAM IDs)
domain_g2_like <- "http://planttfdb.gao-lab.org/multi_align/family/G2-like/domain_aln.fas"
domain_hdzip <- "http://planttfdb.gao-lab.org/multi_align/family/HD-ZIP/domain_aln.fas"
domain_BELL <- "http://planttfdb.gao-lab.org/multi_align/family/WOX/domain_aln.fas"
domain_wus_type <- "http://planttfdb.gao-lab.org/multi_align/family/WOX/domain_aln.fas"
domain_HRT_like <- "http://planttfdb.gao-lab.org/multi_align/family/HRT-like/domain_aln.fas"
domain_NFYB <- "http://planttfdb.gao-lab.org/multi_align/family/NF-YB/domain_aln.fas"
domain_NFYC <- "http://planttfdb.gao-lab.org/multi_align/family/NF-YC/domain_aln.fas"
domain_SAP <- "http://planttfdb.gao-lab.org/multi_align/family/SAP/domain_aln.fas"
domain_trihelix <- "http://planttfdb.gao-lab.org/multi_align/family/Trihelix/domain_aln.fas"
domain_VOZ <- "http://planttfdb.gao-lab.org/multi_align/family/VOZ/domain_aln.fas"

domains <- c(domains, domain_BELL, domain_g2_like, domain_hdzip,
             domain_wus_type, domain_HRT_like, domain_NFYB, domain_NFYC,
             domain_SAP, domain_trihelix, domain_VOZ)

create_profile_HMM <- function(domains = NULL, outdir = NULL) {
  
  f <- lapply(domains, function(x) {
    
    if(startsWith(x, "PF")) {
      url <- paste0('http://pfam.xfam.org/family/', x,
                    '/alignment/seed/format?format=fasta')
    } else {
      url <- x
      x <- gsub("\\/domain.*", "", x)
      x <- gsub(".*\\/", "", x)
    }
    aln_path <- paste0(tempdir(), "/aln.fasta")
    hmm_path <- paste0(tempdir(), "/", x, ".hmm")
    download.file(url = url, destfile = aln_path)
    build <- system2("hmmbuild", args = c(hmm_path, aln_path))
    move <- fs::file_move(hmm_path, paste0(outdir, "/", x, ".hmm"))
    Sys.sleep(2)
    return(move)
  })
  files <- list.files(outdir, pattern = ".hmm")
  status_df <- data.frame(
    domain = domains,
    build_status = str_detect(files, domains)
  )
  return(status_df)
}

create_profile_HMM(domains, 
                   outdir = here::here("inst", "extdata"))
```
