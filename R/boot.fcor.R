#' Fisher-transformed Pearson's correlation: Bootstrap-Based Heterogeneity Test for Between-study Heterogeneity in Random- or Mixed- Effects Model
#'
#' \code{boot.fcor} returns the bootstrap-based tests of the residual heterogeneity in random- or mixed- effects model of Pearson's correlation coefficients transformed with Fisher's r-to-z transformation (z scores).
#'
#' This function returns the test statistics as well as their p-value and significances using (1) Q-test, (2) Bootstrap-Based Heterogeneity Test with Maximum Likelihood (ML), and (3) Bootstrap-Based Heterogeneity Test with Restricted Maximum Likelihood (REML).
#'
#' The results of significances are classified as "sig" or "n.s" based on the cutoff p-value (i.e., alpha level). "sig" means that the between-study heterogeneity is significantly different from zero whereas "n.s" means the between-study heterogeneity is not significantly different from zero. The default alpha level is 0.05.
#'
#' @param n a vector of sample sizes in each of the included studies.
#' @param z a vector of Fisher-transformed Pearson's correlations.
#' @param ttau size of the magnitude to be tested in the alternative hypothesis of the heterogeneity magnitude test. Default to 0.
#' @param model choice of random- or mixed- effects models. Can only be set to \code{"random"}, or \code{"mixed"}.
#' @param mods optional argument to include one or more moderators in the model. \code{mods} is NULL for random-effects model and a dataframe for mixed-effects model. A single moderator can be given as a vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \eqn{k} rows and as many columns as there are moderator variables. See \code{\link[metafor:rma.uni]{rma}} for more details.
#' @param nrep number of replications used in bootstrap simulations. Default to 10^4.
#' @param p_cut cutoff for p-values, which is the alpha level. Default to 0.05.
#' @param boot.include if true, bootstrap simulation results are included in the output (e.g., bootstrap critical values).
#' @param parallel if true, parallel computing using 2 cores will be performed during bootstrapping stage. Otherwise, for loop is used.
#' @param verbose if true, show the progress of boostrapping.
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
#'
#' @references Zuckerman, M. (1994). Behavioral expressions and biosocial bases of sensation seeking. New York, NY: Cambridge University Press.
#' @references Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. URL: http://www.jstatsoft.org/v36/i03/
#'
#' @examples
#' # A meta-analysis of 13 studies studying the correlation
#' # between sensation seeking scores and levels of monoamine oxidase (Zuckerman, 1994).
#'
#' sensation <- boot.heterogeneity:::sensation
#'
#' # n is a list of samples sizes
#' n <- sensation$n
#'
#' # Pearson's correlation
#' r <- sensation$r
#'
#' # Fisher's Transformation
#' z <- 1/2*log((1+r)/(1-r))
#'
#' \dontrun{
#' #' boot.run <- boot.fcor(n, z, model = 'random', p_cut = 0.05)
#' }
#' @export

boot.fcor <- function(n, z, ttau = 0, model = 'random', mods = NULL, nrep = 10^4, p_cut = 0.05, boot.include = FALSE, parallel = FALSE, verbose = FALSE) {

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

  model.f1<-try(metafor::rma(z, vi, mods = mods, tau2=ttau^2, method="ML")) ####NEW!!!!
  model.f2<-try(metafor::rma(z, vi, mods = mods, tau2=ttau^2, method="REML")) ####NEW!!!!
  model.r1<-try(metafor::rma(z, vi, mods = mods, method="ML"))
  model.r2<-try(metafor::rma(z, vi, mods = mods, method="REML"))

  if (sum(!class(model.r2)!="try-error")==0 ){

  bs <- model.r2$beta[,1]
  z_overall <- apply(cbind(1, mods), 1, function(x) sum(bs*x))

  #get predicted effect size for each study #w/ and w/o moderators

  if(verbose){cat("Bootstrapping... \n")}

  if(parallel){
    find.c <- do.call(cbind, pbmcapply::pbmclapply(1:nrep, simulate.z, z_overall=z_overall, ttau=ttau, vi=vi, n=n, mods=mods, mc.cores = 12))
    # parallel::detectCores()-1)
  } else {
    find.c <- matrix(NA, 3, nrep)
    pb <- utils::txtProgressBar(min = 0, max = nrep, style = 3)
    for(i in 1:nrep){
      Sys.sleep(0.01)
      utils::setTxtProgressBar(pb, i)
      find.c[,i] = simulate.z(i, z_overall, ttau, vi, n, mods)
    }
  }
  err.catcher <- sum(colSums(is.na(find.c))!=0)/nrep
  if (err.catcher >0.05){
    warning("Noncovergence rate in simulations is larger than 5%!")
  }

  # We recommend B-REML-LR
  # p-value
  if (model.r1$tau2>=(ttau^2)){
    # One-sided test so the estimated tau has to be larger than lambda,
    # otherwise, we fail to reject the null hypothesis.
    f <- ecdf(na.omit(unlist(find.c)[ c(FALSE,TRUE,FALSE) ]))
    pvalue=1-f( (fitstats(model.r2)-fitstats( model.f2))[1]*2)
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

  #out <- data.frame(stat = c(Q, lllr1, lllr2), p_value = c(Qp, p_lr1, p_lr2), Heterogeneity = c(Qres, res_lr1, res_lr2))
  #(Q, Qp, Qres, lllr1, p_lr1, res_lr1, lllr2, p_lr2, res_lr2)
  #rownames(out) <- c('Qtest', 'boot.ML', 'boot.REML')

  #out <- data.frame(stat = c(Q, Q, lllr1, lllr2), p_value = c(Qp, p_Q, p_lr1, p_lr2), Heterogeneity = c(Qres, res_bootQ, res_lr1, res_lr2))
  #rownames(out) <- c('Qtest', 'boot.Qtest', 'boot.ML', 'boot.REML')

  out <- data.frame(stat = c(Q, lllr2), p_value = c(Qp, p_lr2), Heterogeneity = c(Qres, res_lr2))
  rownames(out) <- c('Qtest', 'boot.REML')

  if(boot.include){
    out <- list(results = out, find.c = find.c, ML.sim = ML.sim, REML.sim = REML.sim, chisq.sim = chisq.sim)
  }

  return(out)
}

