#' Standardized Mean Differences (d): Bootstrap-based Heterogeneity Test for Between-study Heterogeneity in Random- or Mixed- Effects Model
#'
#' \code{boot.d} returns the bootstrap-based tests of the residual heterogeneity in random- or mixed- effects model of standardized mean differences (d).
#'
#' For standardized mean difference, if the biased estimates (i.e., g values) are provided, \code{adjust=TRUE} can be specified to obtain the corresponding unbiased estimates.
#'
#' This function returns the test statistics as well as their p-value and significances using (1) Q-test and (2) Bootstrap-based Heterogeneity Test with Restricted Maximum Likelihood (REML).
#'
#' The results of significances are classified as "sig" or "n.s" based on the cutoff p-value (i.e., alpha level). "sig" means that the between-study heterogeneity is significantly different from zero whereas "n.s" means the between-study heterogeneity is not significantly different from zero. The default alpha level is 0.05.
#'
#' @param n1 A vector of sample sizes from group 1 in each of the included studies.
#' @param n2 A vector of sample sizes from group 2 in each of the included studies.
#' @param est A vector of unbiased estimates of standardized mean differences.
#' @param lambda Size of the magnitude to be tested in the alternative hypothesis of the heterogeneity magnitude test. Default to 0.
#' @param adjust If biased estimates (i.e., g values) are provided, \code{adjust} must be set to \code{TRUE} to compensate for small sample bias. By default, \code{adjust} is set to \code{FALSE}.
#' @param model Choice of random- or mixed- effects models. Can only be set to \code{"random"}, or \code{"mixed"}.
#' @param mods Optional argument to include moderators in the model. \code{mods} is NULL for random-effects model and a dataframe of moderators for mixed-effects model. A single moderator can be given as a vector specifying the values of the moderator. Multiple moderators are specified by giving a matrix with as many columns as there are moderator variables. See \code{\link[metafor:rma.uni]{rma}} for more details.
#' @param nrep Number of replications used in bootstrap simulations. Default to 10^4.
#' @param p_cut Cutoff for p-values, which is the alpha level. Default to 0.05.
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

#' @references Hedges, L. V. (1981). Distribution theory for glass’s estimator of effect size and related estimators. Journal of Educational and Behavioral Statistics, 6(2), 107–128.
#' @references Hedges, L. V., & Olkin, I. (1985). Statistical methods for meta-analysis. San Diego, CA: Academic Press.
#' @references Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. URL: http://www.jstatsoft.org/v36/i03/

#' @examples
#' # Demo 1: A meta-analysis of 18 studies in which the effect of open versus
#' # traditional education on students' self-concept was studied (Hedges & Olkin, 1985).
#'
#' selfconcept <- boot.heterogeneity:::selfconcept
#'
#' # n1 and n2 are lists of samples sizes in two groups
#' n1 <- selfconcept$n1
#' n2 <- selfconcept$n2
#'
#' # g is a list of biased estimates of standardized mean differences in the meta-analytical study
#' g <- selfconcept$g
#' cm <- (1-3/(4*(n1+n2-2)-1)) #correct factor to compensate for small sample bias (Hedges, 1981)
#' d <- cm*g
#'
#' \dontrun{
#' boot.run <- boot.d(n1, n2, est = d, model = 'random', p_cut = 0.05)
#' # is equivalent to:
#' boot.run2 <- boot.d(n1, n2, est = g, model = 'random', adjust = TRUE, p_cut = 0.05)
#' }
#'
#'# Demo 2: A hypothetical meta-analysis of 15 studies with 3 moderators.
#' hypo_moder <- boot.heterogeneity:::hypo_moder
#' \dontrun{
#' boot.run3 <- boot.d(n1 = hypo_moder$n1, n2 = hypo_moder$n2, est = hypo_moder$d, model = 'mixed',
#' mods = cbind(hypo_moder$cov.z1, hypo_moder$cov.z2, hypo_moder$cov.z3), p_cut = 0.05)
#' }
#'
#' @export

boot.d <- function(n1, n2, est, lambda = 0, model = 'random', adjust = FALSE, mods = NULL, nrep = 10^4, p_cut = 0.05, boot.include = FALSE, parallel = FALSE, cores = 4, verbose = FALSE) {

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
  # adjustment for bias
  if(adjust){
    cm <- (1-3/(4*(n1+n2-2)-1)) #correct factor to compensate for small sample bias (Hedges, 1981)
    est <- cm*est
  }
  #########################################################################
  vi<-(n1+n2)/n1/n2+est^2/(2*(n1+n2))

  if(is.null(mods)){
    model.f1<-try(metafor::rma(est, vi, tau2=lambda^2, method="ML"))
    model.f2<-try(metafor::rma(est, vi, tau2=lambda^2, method="REML"))
    model.r1<-try(metafor::rma(est, vi, method="ML"))
    model.r2<-try(metafor::rma(est, vi, method="REML"))
  } else {
    model.f1<-try(metafor::rma(est, vi, mods = mods, tau2=lambda^2, method="ML"))
    model.f2<-try(metafor::rma(est, vi, mods = mods, tau2=lambda^2, method="REML"))
    model.r1<-try(metafor::rma(est, vi, mods = mods, method="ML"))
    model.r2<-try(metafor::rma(est, vi, mods = mods, method="REML"))
  }

  if (sum(!class(model.r2)!="try-error")==0 ){

  bs <- model.r2$beta[,1]
  d_overall <- apply(cbind(1, mods), 1, function(x) sum(bs*x))
  #get predicted effect size for each study #for w/ and w/o moderators

  if(verbose){cat("Bootstrapping... \n")}

  if(parallel){
    find.c <- do.call(cbind, pbmcapply::pbmclapply(1:nrep, simulate.d, d_overall=d_overall, lambda=lambda, vi=vi, n1=n1, n2=n2, mods=mods, mc.cores = cores))
    # parallel::detectCores()-1)
  } else {
    find.c <- matrix(NA, 3, nrep)
    pb <- utils::txtProgressBar(min = 0, max = nrep, style = 3)
    for(i in 1:nrep){
      Sys.sleep(0.01)
      utils::setTxtProgressBar(pb, i)
      find.c[,i] = simulate.d(i, d_overall, lambda, vi, n1, n2, mods)
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
      lllr1<-NA; p_lr1<-NA; res_lr1<-NA; p_lr1.a<-NA; p_Q<-NA;
    }

    if (sum(!class(model.r2)!="try-error" , !class(model.f2)!="try-error")==0){
      lllr2<-(metafor::fitstats(model.r2)-metafor::fitstats(model.f2))[1]*2
      p_lr2<-sum(REML.sim>=lllr2)/length(REML.sim)
      p_lr2.a <-sum(REML.sim>=2.71)/length(REML.sim)
      res_lr2<-ifelse(lllr2>REML.c, 'sig', 'n.s')
    } else {
      lllr2<-NA; p_lr2<-NA; res_lr2<-NA; p_lr2.a<-NA
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

