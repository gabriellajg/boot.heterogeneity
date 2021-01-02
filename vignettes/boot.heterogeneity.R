## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----github, eval = FALSE-----------------------------------------------------
#  # install.packages("devtools")
#  library(devtools)
#  devtools::install_github("gabriellajg/boot.heterogeneity",
#                           force = TRUE,
#                           build_vignettes = TRUE,
#                           dependencies = TRUE)
#  library(boot.heterogeneity)

## ---- eval=FALSE--------------------------------------------------------------
#  library(metafor) # for Q-test
#  library(pbmcapply) # optional - for parallel implementation of bootstrapping
#  library(HSAUR3) # for an example dataset in the tutorial
#  library(knitr) # for knitting the tutorial
#  library(rmarkdown) # for knitting the tutorial

## -----------------------------------------------------------------------------
selfconcept <- boot.heterogeneity:::selfconcept

## -----------------------------------------------------------------------------
head(selfconcept, 3)

## -----------------------------------------------------------------------------
# n1 and n2 are lists of samples sizes in two groups
n1 <- selfconcept$n1
n2 <- selfconcept$n2
# g is a list of effect sizes
g <- selfconcept$g

## -----------------------------------------------------------------------------
cm <- (1-3/(4*(n1+n2-2)-1)) #correct factor to compensate for small sample bias (Hedges, 1981)
d <- cm*g

## ---- eval=FALSE, results = 'hide'--------------------------------------------
#  boot.run <- boot.d(n1, n2, est = d, model = 'random', p_cut = 0.05)

## ---- eval=FALSE, results = 'hide'--------------------------------------------
#  boot.run2 <- boot.d(n1, n2, est = g, model = 'random', adjust = TRUE, p_cut = 0.05)

## ---- eval=FALSE--------------------------------------------------------------
#  boot.run
#  #>                  stat  p_value Heterogeneity
#  #> Qtest       23.391659 0.136929           n.s
#  #> boot.REML    2.037578 0.053100           n.s

## ---- eval=FALSE--------------------------------------------------------------
#  boot.run2
#  #>                  stat  p_value Heterogeneity
#  #> Qtest       23.391659 0.136929           n.s
#  #> boot.REML    2.037578 0.053100           n.s

## -----------------------------------------------------------------------------
hypo_moder <- boot.heterogeneity:::hypo_moder

## -----------------------------------------------------------------------------
head(hypo_moder)

## ---- eval=FALSE, results = 'hide'--------------------------------------------
#  boot.run3 <- boot.d(n1 = hypo_moder$n1,
#                  n2 = hypo_moder$n2,
#                  est = hypo_moder$d,
#                  model = 'mixed',
#                  mods = cbind(hypo_moder$cov.z1, hypo_moder$cov.z2, hypo_moder$cov.z3),
#                  p_cut = 0.05)

## ---- eval=FALSE--------------------------------------------------------------
#  boot.run3
#  #>                  stat    p_value  Heterogeneity
#  #> Qtest       31.849952  0.000806             sig
#  #> boot.REML    9.283428  0.000400             sig

## -----------------------------------------------------------------------------
sensation <- boot.heterogeneity:::sensation

## -----------------------------------------------------------------------------
# n is a list of samples sizes
n <- sensation$n
# Pearson's correlation
r <- sensation$r
# Fisher's Transformation
z <- 1/2*log((1+r)/(1-r))

## ---- eval=FALSE, results = 'hide'--------------------------------------------
#  boot.run.cor <- boot.fcor(n, z, model = 'random', p_cut = 0.05)

## ---- eval=FALSE--------------------------------------------------------------
#  boot.run.cor
#  #>                  stat      p_value    Heterogeneity
#  #> Qtest       29.060970    0.00385868             sig
#  #> boot.REML    6.133111    0.00400882             sig

## ---- eval=FALSE, results = 'hide'--------------------------------------------
#  boot.run.cor2 <- boot.fcor(n, z, lambda=0.08, model = 'random', p_cut = 0.05)

## ---- eval=FALSE--------------------------------------------------------------
#  boot.run.cor2
#  #>                  stat      p_value    Heterogeneity
#  #> boot.REML     2.42325   0.04607372              sig

## -----------------------------------------------------------------------------
library(HSAUR3)
data(smoking)

## -----------------------------------------------------------------------------
# Y1: receive treatment; Y2: stop smoking
n_00 <- smoking$tc - smoking$qc  # not receive treatement yet not stop smoking
n_01 <- smoking$qc # not receive treatement but stop smoking
n_10 <- smoking$tt - smoking$qt # receive treatement but not stop smoking
n_11 <- smoking$qt # receive treatement and stop smoking

## -----------------------------------------------------------------------------
lnOR <- log(n_11*n_00/n_01/n_10)
lnOR

## ---- eval=FALSE, results = 'hide'--------------------------------------------
#  boot.run.lnOR <- boot.lnOR(n_00, n_01, n_10, n_11, model = 'random', p_cut = 0.05)

## ---- eval=FALSE--------------------------------------------------------------
#  boot.run.lnOR
#  #>                  stat    p_value    Heterogeneity
#  #> Qtest       34.873957  0.09050857             n.s
#  #> boot.REML    3.071329  0.03706729             sig

## ---- eval=FALSE, results = 'hide'--------------------------------------------
#  boot.run.lnOR2 <- boot.lnOR(n_00, n_01, n_10, n_11, model = 'random', p_cut = 0.05,
#                              parallel = TRUE, cores = 4)

## ---- eval=FALSE--------------------------------------------------------------
#  boot.run.lnOR2
#  #|=====================================================| 100%, Elapsed 00:41
#  #>                  stat    p_value    Heterogeneity
#  #> Qtest       34.873957  0.09050857             n.s
#  #> boot.REML    3.071329  0.03706729             sig

## -----------------------------------------------------------------------------
sessionInfo()

