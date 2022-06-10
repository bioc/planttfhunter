Data creation
================

# Data in data/

## gsu.rda

This object contains the protein sequences of the algae *Galdieria
sulphuraria*, and it was downloaded from PLAZA Diatoms. Only genes that
have domains of interest (used for classification) were kept for package
size issues.

``` r
gsu <- Biostrings::readAAStringSet(
  "ftp://ftp.psb.ugent.be/pub/plaza/plaza_diatoms_01/Fasta/proteome.selected_transcript.gsu.fasta.gz"
)
names(gsu) <- gsub("\ .*", "", names(gsu))

# The code below must be run after identifying TF domains with annotate_pfam()
data("gsu_annotation")
gsu <- gsu[names(gsu) %in% gsu_annotation$Gene]

outfile <- paste0(tempdir(), "/gsu.fasta")
Biostrings::writeXStringSet(gsu, filepath = outfile)
gsu <- Biostrings::readAAStringSet(outfile)

usethis::use_data(
  gsu, compress = "xz", overwrite = TRUE
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
table$Forbidden <- gsub("\\r\\n\\t\\tor", " or ", table$Forbidden)

classification_scheme <- table
classification_scheme$DBD[3] <- gsub("or", "and", classification_scheme$DBD[3])

usethis::use_data(classification_scheme,
                  compress = "xz", overwrite = TRUE)
```

## score_cutoff.rda

This data frame contains the sequence and domain cutoff for each family.
The values were retrieved from the Supplementary Materials of the
PlantTFDB v2 paper.

``` r
library(here)

# Create a list of PFAM domains
pfam_domains <- readLines(here("inst", "extdata", "PFAM.hmm"))
pfam_domains <- pfam_domains[grep("^NAME", pfam_domains)]
pfam_domains <- gsub(".*PF", "PF", pfam_domains)
    
# Create a list of self-built domains
selfbuilt <- readLines(here("inst", "extdata", "self_built.hmm"))
selfbuilt <- selfbuilt[grep("^NAME", selfbuilt)]
selfbuilt <- gsub("NAME  ", "", selfbuilt)
    

score_cutoff <- data.frame(domain = c(pfam_domains, selfbuilt))

score_cutoff$domaincutoff <- c(
    11, 15, 20, NA, 12.2, 19.7, 21.5, NA, 21, 22, # 1:10
    22.5, NA, 17, 16, NA, NA, NA, 16, NA, 20.7, # 11:20
    NA, 20.8, 22, 23, 30.2, 21.5, 26, 150, 21.3, 22, # 21:30
    23, 22, 20.5, 21, 20.9, NA, 25, NA, 29.5, 22, # 31:40
    21, 20.5, 22, 21.3, 25, NA, 23, NA, 20.2, 20.2, # 41:50
    20.2, 21.3, NA, 20, 24, 27, 24, 50, 150, 100, # 51:60
    NA, NA, NA, 16
)

usethis::use_data(score_cutoff, compress = "xz", overwrite = TRUE)
```

## gsu_families.rda

``` r
gsu_families <- classify_tfs(gsu_annotation)

usethis::use_data(gsu_families, compress = "xz")
```

# Data in inst/extdata

## \*.hmm files (profile HMMs)

The seed alignments for all TF families were downloaded from PFAM and
converted to profile HMMs using the command `hmmbuild` from HMMER.
Domains that were self-built by PlantTFDB were kindly provided by the
database maintainers by email.

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

# Define function to download seed alignments for PFAM domains
download_pfam_seeds <- function(domains = NULL) {
  
  d <- lapply(domains, function(x) {
    url <- paste0('http://pfam.xfam.org/family/', x,
                  '/alignment/seed/format?format=fasta')
    aln_path <- paste0(tempdir(), "/", x, ".fasta")
    if(!file.exists(aln_path)) {
      options(timeout = 5)
      result <- tryCatch({
        downloader::download(url = url, destfile = aln_path)
        res <- TRUE
      },
      error = function(e) return(FALSE),
      warning = function(w) return(FALSE)
      )
      
      if (!result) {
        file.remove(aln_path)
      }
    }
    return(NULL)
  })
  files <- list.files(tempdir(), pattern = ".fasta")
  ds <- sapply(domains, function(x) return(TRUE %in% 
                                             stringr::str_detect(files, x)))
  status <- data.frame(
    domain = domains,
    download_status = ds
  )
  return(status)
}

download_aln <- download_pfam_seeds(domains)


# Define function to build profile HMMs from alignments
build_profile_HMM <- function(domains = NULL) {
  
  f <- lapply(domains, function(x) {
    
    aln_path <- paste0(tempdir(), "/", x, ".fasta")
    hmm_path <- paste0(tempdir(), "/", x, ".hmm")
    build <- system2("hmmbuild", args = c(hmm_path, aln_path))
    sedargs <- paste0("s/aln$/", x, "/g")
    rename_aln <- system2("sed", args = c("-i", sedargs, hmm_path))
    return(NULL)
  })
  files <- list.files(tempdir(), pattern = ".hmm")
  build <- sapply(domains, function(x) return(TRUE %in% 
                                             stringr::str_detect(files, x)))
  status_df <- data.frame(
    domain = domains,
    build_status = build
  )
  return(status_df)
}

build_profile_HMM(domains)

# Combine .hmm files in a single one
combine_domains <- function() {
  files <- list.files(tempdir(), pattern = ".hmm", full.names = TRUE)
  files <- paste0(files, collapse = " ")
  out_hmm <- file.path(tempdir(), "PFAM.hmm")
  system2("cat", args = c(files, " > ", out_hmm))
  return(NULL)
}

combine_domains()
fs::file_move(path = file.path(tempdir(), "PFAM.hmm"),
              new_path = here::here("inst", "extdata", "PFAM.hmm"))
```
