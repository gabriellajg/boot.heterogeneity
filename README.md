
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mc.heterogeneity

<!-- badges: start -->

<!-- badges: end -->

The package mc.heterogeneity implements a Monte Carlo Based
Heterogeneity Test for standardized mean differences (d),
Fisher-transformed Pearson’s correlations (r), and
natural-logarithm-transformed odds ratio (OR) in Meta-Analysis Studies.

Depending on the presence of moderators, this Monte Carlo Based Test can
be implemented in the random or mixed-effects model. This package uses
rma() function from the R package metafor to obtain parameter estimates
and likelihood, so installation of R package metafor is required.

## Installation

You can install the released version of mc.heterogeneity from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("mc.heterogeneity")
```

And the development version from
[GitHub](https://github.com/gabriellajg/mc.heterogeneity) with:

``` r
# install.packages("devtools")
# library(devtools)
devtools::install_github("gabriellajg/mc.heterogeneity")
```

If you have already installed this package locally, you need to override
the previous version with `force` installation:

``` r
devtools::install_github("gabriellajg/mc.heterogeneity", force = TRUE)
```
