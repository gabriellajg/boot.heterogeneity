#' Natural-Logarithm-Transformed Odds Ratio (lnOR): Bootstrap-based Heterogeneity Test for Between-study Heterogeneity in Random- or Mixed- Effects Model
#'
#' \code{boot.lnOR} returns the bootstrap-based tests of the residual heterogeneity in random- or mixed- effects model of natural-logarithm-transformed observed odds ratio (lnOR).
#'
#' For odds ratio, its standard error will be infinite if any one of the four cells in the contingency tables is zero. In this case, Haldane and Anscombe correction is used by adding 0.5 to each cell value (Anscombe, 1956; Haldane, 1940).

#' This function returns the test statistics as well as their p-value and significances using (1) Q-test and (2) Bootstrap-based Heterogeneity Test with Restricted Maximum Likelihood (REML).
#'
#' The results of significances are classified as "sig" or "n.s" based on the cutoff p-value (i.e., alpha level). "sig" means that the between-study heterogeneity is significantly different from zero whereas "n.s" means the between-study heterogeneity is not significantly different from zero. The default alpha level is 0.05.
#'
#' @param n_00 A vector of number of participants who score negatively on both Y1 and Y2 (e.g., mortality cases in the control group).
#' @param n_01 A vector of number of participants who score negatively on Y1 and positively on Y2  (e.g., recovery cases in the control group).
#' @param n_10 A vector of number of participants who score positively on Y1 and negatively on Y2 (e.g., mortality cases in the experimental group).
#' @param n_11 A vector of number of participants who score positively on both Y1 and Y2 (e.g., recovery cases in the experimental group).
#' @param lnOR A vector of natural-logarithm-transformed odds ratio in the included studies, which is calculated as ln(n11*n00/n01/n10)
#' @param lambda Size of the magnitude to be tested in the alternative hypothesis of the heterogeneity magnitude test. Default to 0.
#' @param model Choice of random- or mixed- effects models. Can only be set to \code{"random"}, or \code{"mixed"}.
#' @param mods Optional argument to include moderators in the model. \code{mods} is NULL for random-effects model and a dataframe of moderators for mixed-effects model. A single moderator can be given as a vector specifying the values of the moderator. Multiple moderators are specified by giving a matrix with as many columns as there are moderator variables. See \code{\link[metafor:rma.uni]{rma}} for more details.
#' @param nrep Number of replications used in bootstrap simulations. Default to 10^4.
#' @param p_cut Cutoff for p-value, which is the alpha level. Default to 0.05.
#' @param boot.include If true, bootstrap simulation results are included in the output (e.g., bootstrap critical values).
#' @param parallel If true, parallel computing using 4 cores will be performed during bootstrapping stage. Otherwise, for loop is used.
#' @param cores The number of cores used in the parallel computing. Default to 4.
#' @param verbose If true, show the progress of bootstrapping.
#'
#' @return A dataframe that contains the test statistics ('stat'), p-values ('p_value'), and significances of effect size heterogeneity ("Heterogeneity").
#'
#' @importFrom metafor rma
#' @importFrom metafor fitstats
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats na.omit
#' @importFrom stats ecdf

#' @references Anscombe, F. J. (1956). On estimating binomial response relations. Biometrika, 43(3/4), 461–464.
#' @references Haldane, J. (1940). The mean and variance of| chi 2, when used as a test of homogeneity, when expectations are small. Biometrika, 31(3/4), 346–355.
#' @references Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. URL: http://www.jstatsoft.org/v36/i03/
#' @source Silagy  C, Lancaster  T, Stead  LF, Mant  D, Fowler  G. (2004). Nicotine replacement therapy for smoking cessation. Cochrane Database of Systematic Reviews 2004, Issue 3. Art. No.: CD000146. DOI: 10.1002/14651858.CD000146.pub2.
#'
#' @examples
#' # A meta-analysis consists of 26 studies on nicotine replacement therapy for smoking cessation
#' library(HSAUR3)
#' data(smoking)
#'
#' # Y1: receive treatment; Y2: stop smoking
#' n_00 <- smoking$tc - smoking$qc  # not receive treatement yet not stop smoking
#' n_01 <- smoking$qc # not receive treatement but stop smoking
#' n_10 <- smoking$tt - smoking$qt # receive treatement but not stop smoking
#' n_11 <- smoking$qt # receive treatement and stop smoking
#' lnOR <- log(n_11*n_00/n_01/n_10)
#'
#' \dontrun{
#' boot.run <- boot.lnOR(n_00, n_01, n_10, n_11, model = 'random', p_cut = 0.05)
#' }
#' @export

boot.lnOR <- function(n_00, n_01, n_10, n_11, lambda = 0, model = 'random', mods = NULL, nrep = 10^4, p_cut = 0.05, boot.include = FALSE, parallel = FALSE, cores = 4, verbose = FALSE) {

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
  # zero count correction
  df <- cbind(n_00, n_01, n_10, n_11)
  if(any(df == 0)){
    df <- df + 0.5*(df==0)
    n_00 <- df$n_00
    n_01 <- df$n_01
    n_10 <- df$n_10
    n_11 <- df$n_11
    }
  #########################################################################

  n <- n_00 + n_01 + n_10 + n_11
  lnOR <- log(n_11*n_00/n_01/n_10)
  vi <- 1/n_00+1/n_01+1/n_10+1/n_11

  if(is.null(mods)){
    model.f1<-try(metafor::rma(lnOR, vi, tau2=lambda^2, method="ML")) ####NEW!!!!
    model.f2<-try(metafor::rma(lnOR, vi, tau2=lambda^2, method="REML")) ####NEW!!!!
    model.r1<-try(metafor::rma(lnOR, vi, method="ML"))
    model.r2<-try(metafor::rma(lnOR, vi, method="REML"))
  } else {
    model.f1<-try(metafor::rma(lnOR, vi, mods = mods, tau2=lambda^2, method="ML")) ####NEW!!!!
    model.f2<-try(metafor::rma(lnOR, vi, mods = mods, tau2=lambda^2, method="REML")) ####NEW!!!!
    model.r1<-try(metafor::rma(lnOR, vi, mods = mods, method="ML"))
    model.r2<-try(metafor::rma(lnOR, vi, mods = mods, method="REML"))
  }

  if (sum(!class(model.r2)!="try-error")==0){

  bs <- model.r2$beta[,1]
  lnOR_overall <- apply(cbind(1, mods), 1, function(x) sum(bs*x))
  #get predicted effect size for each study #for w/ and w/o moderators

  if(verbose){cat("Bootstrapping... \n")}

  if(parallel){
    find.c <- do.call(cbind, pbmcapply::pbmclapply(1:nrep, simulate.OR, lnOR_overall = lnOR_overall, lambda=lambda, vi = vi, n = n, n_00_s = n_00, n_01_s = n_01, n_10_s = n_10, n_11_s = n_11, mods = mods, mc.cores = cores))
    # parallel::detectCores()-1)
  } else {
    find.c <- matrix(NA, 3, nrep)
    pb <- utils::txtProgressBar(min = 0, max = nrep, style = 3)
    for(i in 1:nrep){
      Sys.sleep(0.01)
      utils::setTxtProgressBar(pb, i)
      find.c[,i] = simulate.OR(i, lnOR_overall, lambda, vi, n, n_00, n_01, n_10, n_11, mods)
    }
  }
  err.catcher <- sum(colSums(is.na(find.c))!=0)/nrep
  if (err.catcher >0.05){
    warning("Noncovergence rate in simulations is larger than 5%!")
  }

  # We recommend B-REML-LR
  # p-value
  if (model.r1$tau2>=(lambda^2)){
    # One-sided test so the estimated tau has to be larger than lambda,
    # otherwise, we fail to reject the null hypothesis.
    #f <- ecdf(na.omit(unlist(find.c)[ c(FALSE,TRUE,FALSE) ]))
    #pvalue=1-f( (fitstats(model.r2)-fitstats( model.f2))[1]*2)
    # Ge's way to calculate p-value (it's the same)
    ML.sim <- stats::na.omit(unlist(find.c)[ c(TRUE,FALSE,FALSE) ])
  REML.sim <- stats::na.omit(unlist(find.c)[ c(FALSE,TRUE,FALSE) ])
  chisq.sim <- stats::na.omit(unlist(find.c)[ c(FALSE,FALSE,TRUE) ])
  ML.c<-stats::quantile(ML.sim, 0.95)
  REML.c<-stats::quantile(REML.sim, 0.95)
  chisq.c<-stats::quantile(chisq.sim, 0.95)

  if (sum(!class(model.r1)!="try-error" , !class(model.f1)!="try-error")==0){
    lllr1<-(metafor::fitstats(model.r1)-metafor::fitstats(model.f1))[1]*2
    p_lr1<-sum(ML.sim>=lllr1)/length(ML.sim)
    p_lr1.a <-sum(ML.sim>=2.71)/length(ML.sim)
    p_Q <- sum(chisq.sim>=model.f1$QE)/length(chisq.sim)
    res_lr1<-ifelse(lllr1>ML.c, 'sig', 'n.s')
    res_bootQ<-ifelse(model.f1$QE>=chisq.c, 'sig', 'n.s')
  } else {
    lllr1<-NA; p_lr1<-NA; res_lr1<-NA; p_lr1.a<-NA; p_Q<-NA; p_lr2.a<-NA
  }

  if (sum(!class(model.r2)!="try-error" , !class(model.f2)!="try-error")==0){
    lllr2<-(metafor::fitstats(model.r2)-metafor::fitstats(model.f2))[1]*2
    p_lr2<-sum(REML.sim>=lllr2)/length(REML.sim)
    p_lr2.a <-sum(REML.sim>=2.71)/length(REML.sim)
    res_lr2<-ifelse(lllr2>REML.c, 'sig', 'n.s')
  } else {
    lllr2<-NA; p_lr2<-NA; res_lr2<-NA
  }

  Q <- model.f1$QE
  Qp <- model.r2$QEp
  Qres<-ifelse(Qp<= p_cut, 'sig', 'n.s') ### vary the size
  } else {
    pvalue=NA
    warning("pvalue=NA and we fail to reject the null hypothesis.")
    #We don't calculate pvalue in this case and just say that we fail to reject the null hypothesis.
  }
  } else {
    Q<-NA
    Qp<-NA
    Qres<-NA
    p_Q<-NA
    res_bootQ<-NA
    lllr1<-NA
    p_lr1<-NA
    res_lr1<-NA
    lllr2<-NA
    p_lr2<-NA
    res_lr2<-NA
  }

  #out <- data.frame(stat = c(Q, lllr1, lllr2), p_value = c(Qp, p_lr1, p_lr2), Heterogeneity = c(Qres, res_lr1, res_lr2))
  #(Q, Qp, Qres, lllr1, p_lr1, res_lr1, lllr2, p_lr2, res_lr2)
  #rownames(out) <- c('Qtest', 'boot.ML', 'boot.REML')

  #out <- data.frame(stat = c(Q, Q, lllr1, lllr2), p_value = c(Qp, p_Q, p_lr1, p_lr2), Heterogeneity = c(Qres, res_bootQ, res_lr1, res_lr2))
  #rownames(out) <- c('Qtest', 'boot.Qtest', 'boot.ML', 'boot.REML')

  if(lambda==0){
    out <- data.frame(stat = c(Q, lllr2), p_value = c(Qp, p_lr2), Heterogeneity = c(Qres, res_lr2))
    rownames(out) <- c('Qtest', 'boot.REML')
  } else {
    out <- data.frame(stat = c(lllr2), p_value = c(p_lr2), Heterogeneity = c(res_lr2))
    rownames(out) <- c('boot.REML')
  }

  if(boot.include){
    out <- list(results = out, find.c = find.c, ML.sim = ML.sim, REML.sim = REML.sim, chisq.sim = chisq.sim)
  }

  return(out)
}

