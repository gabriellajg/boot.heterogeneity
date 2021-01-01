context("Test p-values")
library(boot.heterogeneity)
library(testthat)

test_that("boot.d", {
  selfconcept <- boot.heterogeneity:::selfconcept
  n1 <- selfconcept$n1
  n2 <- selfconcept$n2
  g <- selfconcept$g
  cm <- (1-3/(4*(n1+n2-2)-1))
  d <- cm*g
  boot.run30 <- boot.d(n1, n2, est = d, nrep = 30)

  expect_error(boot.d(n1, n2, est = d, model = 'cc', nrep = 30), "The meta-analytical model must be either random- or mixed- effects model!")
  expect_equal(round(boot.run30[2,2], 2), 0.1)
})


test_that("boot.fcor", {
  sensation <- boot.heterogeneity:::sensation
  n <- sensation$n
  r <- sensation$r
  z <- 1/2*log((1+r)/(1-r))
  boot.run.cor30 <- boot.fcor(n, z, nrep = 30)

  expect_error(boot.fcor(n, z, model = 'cc', nrep = 30), "The meta-analytical model must be either random- or mixed- effects model!")
  expect_equal(round(boot.run.cor30[2,2], 2), 0)
})


test_that("boot.lnOR", {
  library(HSAUR3)
  data(smoking)
  n_00 <- smoking$tc - smoking$qc  # not receive treatement yet not stop smoking
  n_01 <- smoking$qc # not receive treatement but stop smoking
  n_10 <- smoking$tt - smoking$qt # receive treatement but not stop smoking
  n_11 <- smoking$qt # receive treatement and stop smoking
  lnOR <- log(n_11*n_00/n_01/n_10)
  boot.run.lnOR30 <- boot.lnOR(n_00, n_01, n_10, n_11, nrep = 30)

  expect_error(boot.lnOR(n_00, n_01, n_10, n_11, model = 'cc', nrep = 30), "The meta-analytical model must be either random- or mixed- effects model!")
  expect_equal(round(boot.run.lnOR30[2,2], 2), 0.03)
})


# reference: https://r-pkgs.org/tests.html

