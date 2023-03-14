
<!-- README.md is generated from README.Rmd. Please edit that file -->

# planttfhunter <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/almeidasilvaf/planttfhunter)](https://github.com/almeidasilvaf/planttfhunter/issues)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check-bioc](https://github.com/almeidasilvaf/planttfhunter/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/almeidasilvaf/planttfhunter/actions)
[![Codecov test
coverage](https://codecov.io/gh/almeidasilvaf/planttfhunter/branch/devel/graph/badge.svg)](https://codecov.io/gh/almeidasilvaf/planttfhunter?branch=devel)
<!-- badges: end -->

The goal of **planttfhunter** is to identify plant transcription factors
from protein sequence data and classify them into families and
subfamilies using the classification scheme implemented in
[PlantTFDB](https://doi.org/10.1093/nar/gkz1020).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install **planttfhunter** from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("planttfhunter")
```

And the development version from
[GitHub](https://github.com/almeidasilvaf/planttfhunter) with:

``` r
BiocManager::install("almeidasilvaf/planttfhunter")
```

## Citation

Below is the citation output from using `citation('planttfhunter')` in
R. Please run this yourself to check for any updates on how to cite
**planttfhunter**.

``` r
print(citation('planttfhunter'), bibtex = TRUE)
#> 
#> To cite package 'planttfhunter' in publications use:
#> 
#>   Almeida-Silva F, Van de Peer Y (2022). _planttfhunter: Identification
#>   and classification of plant transcription factors_. R package version
#>   0.99.0, <https://github.com/almeidasilvaf/planttfhunter>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {planttfhunter: Identification and classification of plant transcription factors},
#>     author = {Fabrício Almeida-Silva and Yves {Van de Peer}},
#>     year = {2022},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/almeidasilvaf/planttfhunter},
#>   }
```

Please note that the **planttfhunter** project was only made possible
thanks to many other R and bioinformatics software authors, which are
cited either in the vignettes and/or the paper(s) describing this
package.

## Code of Conduct

Please note that the **planttfhunter** project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.15/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation
  website](http://almeidasilvaf.github.io/planttfhunter) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.15/biocthis)*.
