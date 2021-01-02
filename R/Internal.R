# MC simulation function for Standardized Mean Differences (d)
simulate.d<-function(nrep, d_overall, lambda, vi, n1, n2, mods){
  set.seed(nrep)
  delta.s<-rnorm(length(n1), mean=d_overall, sd=lambda) ####NEW!!!!
  d.s<-stats::rnorm(length(n1),mean=delta.s, sd=sqrt(vi))
  vi.s<-(n1+n2)/n1/n2+d.s^2/(2*(n1+n2))
  model.f1.s<-suppressWarnings(try(metafor::rma(d.s, vi.s, mods = mods, tau2=lambda^2, method="ML"), silent = TRUE))
  model.f2.s<-suppressWarnings(try(metafor::rma(d.s, vi.s, mods = mods, tau2=lambda^2,method="REML"), silent = TRUE))
  model.r1.s<-suppressWarnings(try(metafor::rma(d.s, vi.s, mods = mods, method="ML"), silent = TRUE))
  model.r2.s<-suppressWarnings(try(metafor::rma(d.s, vi.s, mods = mods, method="REML"), silent = TRUE))

  #### NEW below!!!!
  if (sum(!class(model.r1.s)!="try-error" , !class(model.f1.s)!="try-error")==0){
    if (model.r1.s$tau2>=(lambda^2)){
      lllr1.s<-(metafor::fitstats(model.r1.s)-metafor::fitstats( model.f1.s))[1]*2
      chisq<-model.f1.s$QE
    } else {
      lllr1.s<-0}
  } else {
    lllr1.s<-NA
    chisq<-NA
  }

  if (sum(!class(model.r2.s)!="try-error", !class(model.f2.s)!="try-error")==0){
    if (model.r2.s$tau2>=(lambda^2)){
      lllr2.s<-(metafor::fitstats(model.r2.s)-metafor::fitstats(model.f2.s))[1]*2
      chisq<-model.f2.s$QE
    } else {
      lllr2.s<-0
      chisq<-length(n1)-1}
  } else {
    lllr2.s<-NA
    chisq<-NA}

  #Sys.sleep(runif(1))
  return(c(lllr1.s, lllr2.s, chisq))
}

# MC simulation function for Fisher's Transformed z Scores (r, z)
simulate.z<-function(nrep, z_overall, lambda, vi, n, mods){
  set.seed(nrep)
  delta.s<-rnorm(length(n), mean=z_overall, sd=lambda) ####NEW!!!!
  z.s<-stats::rnorm(length(n), mean = delta.s, sd = sqrt(vi))
  vi.s<- vi
  model.f1.s<-suppressWarnings(try(metafor::rma(z.s, vi.s, mods = mods, tau2=lambda^2, method="ML"), silent = TRUE))
  model.f2.s<-suppressWarnings(try(metafor::rma(z.s, vi.s, mods = mods, tau2=lambda^2,method="REML"), silent = TRUE))
  model.r1.s<-suppressWarnings(try(metafor::rma(z.s, vi.s, mods = mods, method="ML"), silent = TRUE))
  model.r2.s<-suppressWarnings(try(metafor::rma(z.s, vi.s, mods = mods, method="REML"), silent = TRUE))

  #### NEW below!!!!
  if (sum(!class(model.r1.s)!="try-error" , !class(model.f1.s)!="try-error")==0){
    if (model.r1.s$tau2>=(lambda^2)){
      lllr1.s<-(metafor::fitstats(model.r1.s)-metafor::fitstats( model.f1.s))[1]*2
      chisq<-model.f1.s$QE
    } else {
      lllr1.s<-0}
  } else {
    lllr1.s<-NA
    chisq<-NA
  }

  if (sum(!class(model.r2.s)!="try-error", !class(model.f2.s)!="try-error")==0){
    if (model.r2.s$tau2>=(lambda^2)){
      lllr2.s<-(metafor::fitstats(model.r2.s)-metafor::fitstats(model.f2.s))[1]*2
      chisq<-model.f2.s$QE
    } else {
      lllr2.s<-0
      chisq<-length(n)-1}
  } else {
    lllr2.s<-NA
    chisq<-NA}

  return(c(lllr1.s,lllr2.s, chisq))
}

# MC simulation function for log odds ratio

#n4=dat$qt;   nt=dat$tt; n2=dat$qc; nc=dat$tc;n3=nt-n4; n1=nc-n2
#n_11                    n_01                 n_10      n_00

simulate.OR<-function(nrep, lnOR_overall, lambda, vi, n, n_00_s, n_01_s, n_10_s, n_11_s, mods){
  set.seed(nrep)
  delta.s <- stats::rnorm(length(n), mean=lnOR_overall, sd=lambda)
  lnOR.s <- stats::rnorm(length(n), mean=delta.s, sd=sqrt(vi))
  # n_00_s <- n_00
  # n_01_s <- n_01
  # n_10_s <- n_00_s*(n-n_00_s-n_01_s)/(n_00_s + n_01_s*exp(lnOR.s))
  # n_11_s <- n - n_00_s - n_01_s - n_10_s

  index<-sample(1:4,1,prob=c(0.25,0.25,0.25,0.25))

  if (index==1){
    n_00_s<-exp(lnOR.s)*n_10_s*n_01_s/n_11_s
  }

  if (index==2){
    n_01_s<-n_11_s*n_00_s/exp(lnOR.s)/n_10_s
  }

  if (index==3){
    n_10_s<-n_11_s*n_00_s/exp(lnOR.s)/n_01_s
  }

  if (index==4){
    n_11_s<-exp(lnOR.s)*n_10_s*n_01_s/n_00_s
  }

  #########################################################################
  vi.s <- 1/n_00_s+1/n_01_s+1/n_10_s+1/n_11_s
  lnOR.s <- log(n_11_s*n_00_s/n_01_s/n_10_s) ###new!!

  # zero count correction
  df <- cbind(n_00_s, n_01_s, n_10_s, n_11_s)
  if(any(df == 0)){
    df <- df + 0.5*(df==0)
    n_00_s <- df$n_00_s
    n_01_s <- df$n_01_s
    n_10_s <- df$n_10_s
    n_11_s <- df$n_11_s
  }
  #########################################################################
  model.f1.s<-suppressWarnings(try(metafor::rma(lnOR.s, vi.s, mods = mods, tau2=lambda^2, method="ML"), silent = TRUE))
  model.f2.s<-suppressWarnings(try(metafor::rma(lnOR.s, vi.s, mods = mods, tau2=lambda^2,method="REML"), silent = TRUE))
  model.r1.s<-suppressWarnings(try(metafor::rma(lnOR.s, vi.s, mods = mods, method="ML"), silent = TRUE))
  model.r2.s<-suppressWarnings(try(metafor::rma(lnOR.s, vi.s, mods = mods, method="REML"), silent = TRUE))

  #### NEW below!!!!
  if (sum(!class(model.r1.s)!="try-error" , !class(model.f1.s)!="try-error")==0){
    if (model.r1.s$tau2>=(lambda^2)){
      lllr1.s<-(metafor::fitstats(model.r1.s)-metafor::fitstats( model.f1.s))[1]*2
      chisq<-model.f1.s$QE
    } else {
      lllr1.s<-0}
  } else {
    lllr1.s<-NA
    chisq<-NA
  }

  if (sum(!class(model.r2.s)!="try-error", !class(model.f2.s)!="try-error")==0){
    if (model.r2.s$tau2>=(lambda^2)){
      lllr2.s<-(metafor::fitstats(model.r2.s)-metafor::fitstats(model.f2.s))[1]*2
      chisq<-model.f2.s$QE
    } else {
      lllr2.s<-0
      chisq<-length(n)-1}
  } else {
    lllr2.s<-NA
    chisq<-NA}

  return(c(lllr1.s,lllr2.s, chisq))
}
