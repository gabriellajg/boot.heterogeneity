
<!-- README.md is generated from README.Rmd. Please edit that file -->

# boot.heterogeneity

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/gabriellajg/boot.heterogeneity.svg?branch=master)](https://travis-ci.org/gabriellajg/boot.heterogeneity) [![](http://cranlogs.r-pkg.org/badges/grand-total/boot.heterogeneity)](https://CRAN.R-project.org/package=boot.heterogeneity)
<!-- badges: end -->

The package boot.heterogeneity implements a Bootstrap-Based
Heterogeneity Test for standardized mean differences (d),
Fisher-transformed Pearson’s correlations (r), and
natural-logarithm-transformed odds ratio (OR) in Meta-Analysis Studies.

Depending on the presence of moderators, this Bootstrap-Based Test can
be implemented in the random or mixed-effects model. This package uses
rma() function from the R package metafor to obtain parameter estimates
and likelihood, so installation of R package metafor is required.

## Installation

You can install the released version of boot.heterogeneity from
[CRAN](https://cran.r-project.org/package=mc.heterogeneity) with:

``` r
install.packages("boot.heterogeneity")
```

And the development version from
[GitHub](https://github.com/gabriellajg/boot.heterogeneity) with:

``` r
# install.packages("devtools")
# library(devtools)
devtools::install_github("gabriellajg/boot.heterogeneity")
```

If you have already installed this package locally, you need to override
the previous version with `force` installation:

``` r
devtools::install_github("gabriellajg/boot.heterogeneity", force = TRUE)
```

Updates to this package will always be pushed to Github first. If you
encounter any errors, please re-install the package from github using
`force` installation and try again.
