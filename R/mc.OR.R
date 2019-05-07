#' Monte Carlo Based Heterogeneity Test for Residual Variances in Random- or Mixed- Effects Model of Log Odds Ratio (OR)
#'
#' \code{mc.d} returns the Monte Carlo based tests of the residual heterogeneity in random- or mixed- effects model of natural-logarithm-transformed observed odds ratio (OR).
#'
#' This function returns the test statistics as well as their significances using (1) Q-test, (2) Monte Carlo Based Heterogeneity Test with Maximum Likelihood (ML), and (3) Monte Carlo Based Heterogeneity Test with Restricted Maximum Likelihood (REML).
#'
#' The results of significances are classified as "sig" or "n.s" based on the cutoff p-value. "sig" means that the residual variance is significantly different from zero wheras "n.s" means the residual variance is not significantly different from zero.
#'
#' @param n a vector of sample sizes in each of the included studies.
#' @param LOR a vector of log odds ratios in each of the included studies.
#' @param p_00 a vector of probabilities of negatives on both Y1 and Y2.
#' @param p_01 a vector of probabilities of negative on Y1 and positive on Y2.
#' @param model choice of random- or mixed- effects models. Can only be set to \code{"random"}, or \code{"mixed"}.
#' @param mods optional argument to include one or more moderators in the model. \code{mods} is NULL for random-effects model and a dataframe for mixed-effects model. A single moderator can be given as a vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \eqn{k} rows and as many columns as there are moderator variables. See \code{\link[metafor]{rma}} for details.
#' @param nrep number of replications used in Monte Carlo Simulations. Default to 10^4.
#' @param p_cut cutoff for p-values. Default to 0.05.
#'
#' @examples
#' library(HSAUR2)
#' data(smoking)
#' # Total number of subjecs in each study
#' n <- smoking$tt + smoking$tc
#' # Y1: receive treatment; Y2: stop smoking
#' p_00 <- (smoking$tc - smoking$qc)/n  # not receive treatement yet not stop smoking
#' p_01 <- smoking$qc/n # not receive treatement but stop smoking
#' p_10 <- (smoking$tt - smoking$qt)/n # receive treatement but not stop smoking
#' p_11 <- smoking$qt/n # receive treatement and stop smoking
#' LOR <- log(p_11*p_00/p_01/p_10)
#' \dontrun{
#' mc.run <- mc.OR(n, LOR, p_00, p_01, model = 'random')
#' }
#' @export

mc.OR <- function(n, LOR, p_00, p_01, model = 'random', mods = NULL, nrep = 10^4, p_cut = 0.05) {

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
  p_10 <- (1-p_00-p_01)/(1 + p_01*exp(LOR)/p_00)
  p_11 <- 1-p_00-p_01-p_10
  n_00 <- n*p_00
  n_01 <- n*p_01
  n_10 <- n*p_10
  n_11 <- n*p_11
  vi <- 1/n_00+1/n_01+1/n_10+1/n_11

  model.f1<-try(metafor::rma(LOR, vi, mods = mods, tau2=0, method="ML"))
  model.f2<-try(metafor::rma(LOR, vi, mods = mods, tau2=0, method="REML"))
  model.r1<-try(metafor::rma(LOR, vi, mods = mods, method="ML"))
  model.r2<-try(metafor::rma(LOR, vi, mods = mods, method="REML"))


  #if (class(model.r2)!="try-error" ){
  if (sum(!class(model.r2)!="try-error")==0){

  bs <- model.r2$beta[,1]
  lor_overall <- apply(cbind(1, mods), 1, function(x) sum(bs*x))
  #get predicted effect size for each study #for w/ and w/o moderators

  find.c <- matrix(NA, 2, nrep)
  pb <- txtProgressBar(min = 0, max = nrep, style = 3)
  for(i in 1:nrep){
    Sys.sleep(0.01)
    setTxtProgressBar(pb, i)
    find.c[,i] = simulate.OR(i, lor_overall, vi, n, p_00, p_01, mods)
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
