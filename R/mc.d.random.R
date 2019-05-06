#' Monte Carlo Based Heterogeneity Test for Residual Variances in Random- or Mixed- Effects Model of Standardized Mean Differences (d)
#'
#' \code{mc.d} returns the Monte Carlo based tests of the residual heterogeneity in random- or mixed- effects model of standardized mean differences (d).
#'
#' This function returns the test statistics as well as their significances using (1) Q-test, (2) Monte Carlo Based Heterogeneity Test with Maximum Likelihood (ML), and (3) Monte Carlo Based Heterogeneity Test with Restricted Maximum Likelihood (REML).
#'
#' The results of significances are classified as "sig" or "n.s" based on the cutoff p-value. "sig" means that the residual variance is significantly different from zero wheras "n.s" means the residual variance is not significantly different from zero.
#'
#' @param n1 a vector of sample sizes from group 1 in each study.
#' @param n2 a vector of sample sizes from group 2 in each study.
#' @param d a vector of bias-corrected estimate of standardized mean differences (often reported as hedge's g)
#' @param nrep number of replications used in Monte Carlo Simulations. Default to 10^4.
#' @param p_cut cutoff for p-values. Default to 0.05.
#' @param model choice of random- or mixed- effects models. Can only be set to \code{"random"}, or \code{"mixed"}.
#' @param mods optional argument to include one or more moderators in the model. \code{mods} is NULL for random-effects model and a dataframe for mixed-effects model. A single moderator can be given as a vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \eqn{k} rows and as many columns as there are moderator variables. See \code{\link[metafor]{rma}} for details.
#'

#' @examples
#' #** n1 and n2 are lists of samples sizes in two groups
#' n1 <- c(100, 131, 40, 40, 97, 28, 60, 72, 87, 80, 79, 70, 36, 9, 14, 21, 133, 83)
#' n2 <- c(180, 138, 40, 40, 47, 61, 55, 102, 45, 49, 55, 109, 93, 18, 16, 22, 124, 45)
#' g is a collection of standardized mean differences in the meta-analytical study
#' g <- c(0.100, -0.162, -0.090, -0.049, -0.046, -0.010, -0.431, -0.261, 0.134, 0.019, 0.175, 0.056, 0.045, 0.103, 0.121, -0.482, 0.290, 0.342)
#' cm <- (1-3/(4*(n1+n2-2)-1)) #correct factor to compensate for small sample bias (Hedges & Olkin, 1985)
#' d <- cm*g
#' mc.run <- mc.d(n1, n2, d, model = 'random')
#' mc.run2 <- mc.d(n1, n2, d, model = 'mixed', mods = tt, nrep = 1000)
#' @export

mc.d <- function(n1, n2, d, model = 'random', mods = NULL, nrep = 10^4, p_cut = 0.05) {

  #########################################################################
  if (!model %in% c('random', 'mixed')){
    stop("The meta-analytical model must be either random- or mixed- effects model!")
  }
  if (model == 'random' & !is.null(mods)){
    stop("No moderators should be included for random-effects model!")
  }
  if (model == 'mixed' & is.null(mods)){
    stop("Moderators need be included for mixed-effects model!")
  }

  #########################################################################
  vi<-(n1+n2)/n1/n2+d^2/(2*(n1+n2))

  model.f1<-try(metafor::rma(d, vi, mods = mods, tau2=0, method="ML"))
  model.f2<-try(metafor::rma(d, vi, mods = mods, tau2=0, method="REML"))
  model.r1<-try(metafor::rma(d, vi, mods = mods, method="ML"))
  model.r2<-try(metafor::rma(d, vi, mods = mods, method="REML"))


  #if (class(model.r2)!="try-error" ){
  if (sum(!class(model.r2)!="try-error")==0 ){

  bs <- model.r2$beta[,1]
  d_overall <- sum(bs*colMeans(cbind(1, mods))) #for w/ and w/o moderators

  find.c <- matrix(NA, 2, nrep)
  pb <- txtProgressBar(min = 0, max = nrep, style = 3)
  for(i in 1:nrep){
    Sys.sleep(0.01)
    setTxtProgressBar(pb, i)
    find.c[,i] = simulate.c(i, d_overall, vi, n1, n2, mods)
    }

  ML.c<-quantile(na.omit(unlist(find.c)[ c(TRUE,FALSE) ]),0.95)
  REML.c<-quantile(na.omit(unlist(find.c)[ c(FALSE,TRUE) ]),0.95)

  if (sum(!class(model.r1)!="try-error" , !class(model.f1)!="try-error")==0){
      lllr1<-(metafor::fitstats(model.r1)-metafor::fitstats(model.f1))[1]*2
      #p_lr1<-lllr1>ML.c
      p_lr1<-ifelse(lllr1>ML.c, 'sig', 'n.s')
  } else {
    lllr1<-NA; p_lr1<-NA
    }

  if (sum(!class(model.r2)!="try-error" , !class(model.f2)!="try-error")==0){
      lllr2<-(metafor::fitstats(model.r2)-metafor::fitstats( model.f2))[1]*2
      #p_lr2<-lllr2>REML.c
      p_lr2<-ifelse(lllr2>REML.c, 'sig', 'n.s')
  } else {
    lllr2<-NA; p_lr2<-NA
    }

  Q <- model.f1$QE
  #p_Q<-model.f1$QEp< p_cut ### vary the size
  p_Q<-ifelse(model.f1$QEp< p_cut, 'sig', 'n.s') ### vary the size
  } else {
    Q<-NA
    lllr1<-NA
    lllr2<-NA
    p_lr1<-NA
    p_lr2<-NA
    p_Q<-NA
  }

  out <- data.frame(Q, p_Q, lllr1, p_lr1, lllr2, p_lr2)
  colnames(out) <- c('QE', 'QEp', 'mc.ML', 'MLp', 'mc.REML', 'REMLp')
  rownames(out) <- NULL
  return(out)
}
