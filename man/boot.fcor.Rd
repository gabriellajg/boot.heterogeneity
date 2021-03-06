% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot.fcor.R
\name{boot.fcor}
\alias{boot.fcor}
\title{Fisher-transformed Pearson's correlation: Bootstrap-based Heterogeneity Test for Between-study Heterogeneity in Random- or Mixed- Effects Model}
\usage{
boot.fcor(
  n,
  z,
  lambda = 0,
  model = "random",
  mods = NULL,
  nrep = 10^4,
  p_cut = 0.05,
  boot.include = FALSE,
  parallel = FALSE,
  cores = 4,
  verbose = FALSE
)
}
\arguments{
\item{n}{A vector of sample sizes in each of the included studies.}

\item{z}{A vector of Fisher-transformed Pearson's correlations.}

\item{lambda}{Size of the magnitude to be tested in the alternative hypothesis of the heterogeneity magnitude test. Default to 0.}

\item{model}{Choice of random- or mixed- effects models. Can only be set to \code{"random"}, or \code{"mixed"}.}

\item{mods}{Optional argument to include moderators in the model. \code{mods} is NULL for random-effects model and a dataframe of moderators for mixed-effects model. A single moderator can be given as a vector specifying the values of the moderator. Multiple moderators are specified by giving a matrix with as many columns as there are moderator variables. See \code{\link[metafor:rma.uni]{rma}} for more details.}

\item{nrep}{Number of replications used in bootstrap simulations. Default to 10^4.}

\item{p_cut}{Cutoff for p-value, which is the alpha level. Default to 0.05.}

\item{boot.include}{If true, bootstrap simulation results are included in the output.}

\item{parallel}{If true, parallel computing using 4 cores will be performed during bootstrapping stage. Otherwise, for loop is used.}

\item{cores}{The number of cores used in the parallel computing. Default to 4.}

\item{verbose}{If true, show the progress of bootstrapping.}
}
\value{
A dataframe that contains the test statistics ('stat'), p-values ('p_value'), and significances of effect size heterogeneity ("Heterogeneity").
}
\description{
\code{boot.fcor} returns the bootstrap-based tests of the residual heterogeneity in random- or mixed- effects model of Pearson's correlation coefficients transformed with Fisher's r-to-z transformation (z scores).
}
\details{
This function returns the test statistics as well as their p-value and significances using (1) Q-test and (2) Bootstrap-based Heterogeneity Test with Restricted Maximum Likelihood (REML).

The results of significances are classified as "sig" or "n.s" based on the cutoff p-value (i.e., alpha level). "sig" means that the between-study heterogeneity is significantly different from zero whereas "n.s" means the between-study heterogeneity is not significantly different from zero. The default alpha level is 0.05.
}
\examples{
# A meta-analysis of 13 studies studying the correlation
# between sensation-seeking scores and levels of monoamine oxidase (Zuckerman, 1994).

sensation <- boot.heterogeneity:::sensation

# n is a list of samples sizes
n <- sensation$n

# Pearson's correlation
r <- sensation$r

# Fisher's Transformation
z <- 1/2*log((1+r)/(1-r))

\dontrun{
#' boot.run <- boot.fcor(n, z, model = 'random', p_cut = 0.05)
}
}
\references{
Zuckerman, M. (1994). Behavioral expressions and biosocial bases of sensation-seeking. New York, NY: Cambridge University Press.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. URL: http://www.jstatsoft.org/v36/i03/
}
