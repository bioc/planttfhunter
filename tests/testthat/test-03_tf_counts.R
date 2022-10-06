
#----Load data------------------------------------------------------------------
set.seed(123)
data(gsu)
data(tf_counts)

# Pick random subsets of 100 genes to simulate other species
proteomes <- list(
    Gsu1 = gsu[sample(names(gsu), 50, replace = FALSE)],
    Gsu2 = gsu[sample(names(gsu), 50, replace = FALSE)],
    Gsu3 = gsu[sample(names(gsu), 50, replace = FALSE)],
    Gsu4 = gsu[sample(names(gsu), 50, replace = FALSE)]
)

# Create species metadata
species_metadata <- data.frame(
    row.names = names(proteomes),
    Division = "Rhodophyta",
    Origin = c("US", "Belgium", "China", "Brazil")
)


#----Start tests----------------------------------------------------------------
test_that("get_tf_counts() returns a SummarizedExperiment object", {
    
    se <- tf_counts
    if(hmmer_is_installed()) {
        se <- get_tf_counts(proteomes, species_metadata)
    }

    expect_error(
        get_tf_counts(proteomes = data.frame(), species_metadata)
    )
    
    expect_error(
        get_tf_counts(proteomes, species_metadata[1:2, ])
    )
    
    expect_true("SummarizedExperiment" %in% class(se))
})
