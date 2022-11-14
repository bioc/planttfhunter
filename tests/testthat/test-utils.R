
#---Load and create data for testing--------------------------------------------
## Example output of hmmsearch (.txt format)
example_output <- c(
    "#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord", 
    "# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target", 
    "#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------", 
    "gsu06140.1           -            340 G2-like              -             56   5.5e-09   24.7   0.0   1   1   2.1e-09   1.1e-08   23.8   0.0     1    53    95   141    95   144 0.89 -", 
    "gsu06160.1           -            120 NF-YC                -            102   1.1e-07   20.5   0.0   1   1   2.7e-08   1.3e-07   20.2   0.0    19    90    18    89     7    98 0.90 -"
)
file <- file.path(tempdir(), "example_file_unit_test.txt")
example_file <- writeLines(
    example_output, file
)

#----Start tests----------------------------------------------------------------
test_that("read_hmmsearch() reads and parses hmmsearch output", {
    output <- read_hmmsearch(file)
    
    expect_equal(ncol(output), 23)
    expect_equal(nrow(output), 2)
    expect_equal(class(output), "data.frame")
})

test_that("filter_hmmer() filters hmmsearch output", {
    output <- read_hmmsearch(file)
    foutput <- filter_hmmer(output)
    
    expect_equal(ncol(foutput), 24)
    expect_equal(nrow(foutput), 1)
    expect_equal(class(foutput), "data.frame")
    expect_equal(
        names(foutput), 
        c(
            "query_name", "domain_name", "domain_accession", "domain_len", 
            "query_accession", "qlen", "sequence_evalue", "sequence_score", 
            "sequence_bias", "domain_N", "domain_of", "domain_cevalue", 
            "domain_ievalue", "domain_score", "domain_bias", "hmm_from", 
            "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc", 
            "description", "domaincutoff"
        )
    )
})

test_that("is_valid_cmd() works", {
    v <- is_valid_cmd(cmd = "hmmsearch", args = "-h")
    expect_equal(class(v), "logical")
    expect_equal(length(v), 1)
})

test_that("hmmer_is_installed() works", {
    h <- hmmer_is_installed()
    expect_equal(class(h), "logical")
    expect_equal(length(h), 1)
})
