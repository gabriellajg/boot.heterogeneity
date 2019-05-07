#' Monte Carlo Based Heterogeneity Test for Residual Variances in Random- or Mixed- Effects Model of Fisher's Transformed z Scores (z)
#'
#' \code{mc.d} returns the Monte Carlo based tests of the residual heterogeneity in random- or mixed- effects model of correlation coefficients transformed with Fisher's r-to-z transformation (z scores).
#'
#' This function returns the test statistics as well as their significances using (1) Q-test, (2) Monte Carlo Based Heterogeneity Test with Maximum Likelihood (ML), and (3) Monte Carlo Based Heterogeneity Test with Restricted Maximum Likelihood (REML).
#'
#' The results of significances are classified as "sig" or "n.s" based on the cutoff p-value. "sig" means that the residual variance is significantly different from zero wheras "n.s" means the residual variance is not significantly different from zero.
#'
#' @param n a vector of sample sizes in each of the included studies.
#' @param z a vector of Fisher's transformed z scores (z).
#' @param model choice of random- or mixed- effects models. Can only be set to \code{"random"}, or \code{"mixed"}.
#' @param mods optional argument to include one or more moderators in the model. \code{mods} is NULL for random-effects model and a dataframe for mixed-effects model. A single moderator can be given as a vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \eqn{k} rows and as many columns as there are moderator variables. See \code{\link[metafor]{rma}} for details.
#' @param nrep number of replications used in Monte Carlo Simulations. Default to 10^4.
#' @param p_cut cutoff for p-values. Default to 0.05.
#'

#' @examples
#' # n is a list of samples sizes in the meta-analytical study.
#' n <- c(65, 30, 93, 36, 57, 30, 40, 13, 44, 58, 125, 10, 13)
#' # Pearson's Correction
#' r <- c(0.17, -0.45, -0.47, -0.13, -0.24, -0.15, -0.25, -0.66, -0.25, -0.23, -0.18, 0.18, -0.74)
#' # Fisher's Transformation
#' z <- 1/2*log((1+r)/(1-r))
#' \dontrun{
#' mc.run <- mc.z(n, z, model = 'random')
#' }
#' @export

mc.z <- function(n, z, model = 'random', mods = NULL, nrep = 10^4, p_cut = 0.05) {

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
  vi<-1/(n-3)

  model.f1<-try(metafor::rma(d, vi, mods = mods, tau2=0, method="ML"))
  model.f2<-try(metafor::rma(d, vi, mods = mods, tau2=0, method="REML"))
  model.r1<-try(metafor::rma(d, vi, mods = mods, method="ML"))
  model.r2<-try(metafor::rma(d, vi, mods = mods, method="REML"))


  #if (class(model.r2)!="try-error" ){
  if (sum(!class(model.r2)!="try-error")==0 ){

  bs <- model.r2$beta[,1]
  z_overall <- apply(cbind(1, mods), 1, function(x) sum(bs*x))
  #get predicted effect size for each study #for w/ and w/o moderators

  find.c <- matrix(NA, 2, nrep)
  pb <- txtProgressBar(min = 0, max = nrep, style = 3)
  for(i in 1:nrep){
    Sys.sleep(0.01)
    setTxtProgressBar(pb, i)
    find.c[,i] = simulate.z(i, z_overall, vi, n, mods)
    }
  err.catcher <- sum(colSums(is.na(find.c))!=0)/nrep
  if (err.catcher >0.05){
    warning("Noncovergence rate in simulations is larger than 5%!")
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
