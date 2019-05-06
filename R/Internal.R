simulate.c<-function(nrep, d_overall, vi, n1, n2, mods){
  set.seed(nrep)
  d.s<-rnorm(length(n1),mean=d_overall, sd=sqrt(vi))
#d.s<-rnorm(length(n1), mean = d_overall[1]+ d_overall[2]*cov.z1+ d_overall[3]*cov.z2+ d_overall[4]*cov.z3, sd=sqrt(vi))
  vi.s<-(n1+n2)/n1/n2+d.s^2/(2*(n1+n2))
  model.f1.s<-try(metafor::rma(d.s, vi.s, mods = mods, tau2=0, method="ML"))
  model.f2.s<-try(metafor::rma(d.s, vi.s, mods = mods, tau2=0,method="REML"))
  model.r1.s<-try(metafor::rma(d.s, vi.s, mods = mods, method="ML"))
  model.r2.s<-try(metafor::rma(d.s, vi.s, mods = mods, method="REML"))

  #if (class(model.r1.s)!="try-error" & class(model.f1.s)!="try-error"){
  if (sum(!class(model.r1.s)!="try-error" , !class(model.f1.s)!="try-error")==0){
    lllr1.s<-(metafor::fitstats(model.r1.s)-metafor::fitstats( model.f1.s))[1]*2} else {
      lllr1.s<-NA}
  #if (class(model.r2.s)!="try-error" & class(model.f2.s)!="try-error"){
  if (sum(!class(model.r2.s)!="try-error" , !class(model.f2.s)!="try-error")==0){
    lllr2.s<-(metafor::fitstats(model.r2.s)-metafor::fitstats( model.f2.s))[1]*2} else {
      lllr2.s<-NA}

  return(c(lllr1.s,lllr2.s))
}

#  mc.d.random: find.c<- mclapply(1:10^4, simulate.c, mc.cores=ncores)
