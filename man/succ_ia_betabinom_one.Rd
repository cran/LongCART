\name{succ_ia_betabinom_one}
\alias{succ_ia_betabinom_one}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determines predictive power of success based on interim results and beta prior for one-sample binary data
}
\description{
This function can be used to determine predictive power for trial success and clinical success based on the interim results and beta prior distribution for test of population proportion. 
}
\usage{
succ_ia_betabinom_one(N, n, x, 
                      null.value = 0, alternative = "greater", 
                      test="z", correct=TRUE,
                      succ.crit = "trial", Z.crit.final = 1.96, 
                      alpha.final = 0.025, clin.succ.threshold = NULL, 
                      a = 1, b = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Sample size at final analysis. Cannot be missing.}
  \item{n}{Sample size at interim analysis. Cannot be missing.}
  \item{x}{Number of observed response at interim analysis. Cannot be missing.}
  \item{null.value}{The specified value under null hypothesis. Default is 0.}
  \item{alternative}{Direction of alternate hypothesis. Can be "greater" or "less".}
  \item{test}{Statistical test. Default is \code{"z"} for Z test. For Exact binomial test, specify \code{"exact"}.}
  \item{correct}{A logical indicating whether Yates' continuity correction should be applied where possible. Applies to approximate Z-test only.}
  \item{succ.crit}{Specify "trial" for trial success (i.e., null hypothesis is rejected at final analysis) or "clinical" for clinical success (i.e., estimated value at the final analysis is greater than clinically meaningful value as specified under \code{clin.succ.threshold}.)}
  \item{Z.crit.final}{The rejection boundary at final analysis in Z-value scale. Either \code{alpha.final} or \code{Z.crit.final} must be specified when determining trial success.}
  \item{alpha.final}{The rejection boundary at final analysis in alpha (1-sided) scale (e.g., 0.025). Either \code{alpha.final} or \code{Z.crit.final} must be specified when determining trial success.}
  \item{clin.succ.threshold}{Clinically meaningful value. Required when \code{succ.crit="clinical"}.}
  \item{a}{Value of \code{a} corresponding to Beta(a, b) prior for proportion.}
  \item{b}{Value of \code{b} corresponding to Beta(a, b) prior for proportion.}
}
\details{
This function can be used to determine Predictive power or predictive probability of success (PPoS) based on the interim results for testing of population proportion. The calculation of PoS is carried out assuming beta prior distributions for proportion. This function can be used to determine clinical success (\code{succ.crit="clinical"}) and trial success (\code{succ.crit="trial"}). For clinical success, \code{clin.succ.threshold} must be specified. For trial success, \code{Z.crit.final or alpha.final} must be specified.
}
\author{Madan Gopal Kundu <madan_g.kundu@yahoo.com>}
\references{
Kundu, M. G., Samanta, S., and Mondal, S. (2021). An introduction to the determination of the probability of a successful trial: Frequentist and Bayesian approaches. arXiv preprint arXiv:2102.13550.
}
\seealso{\code{\link{succ_ia_betabinom_two}}, \code{\link{succ_ia}}, \code{\link{PoS}}   }

\examples{
succ_ia_betabinom_one( N=40, n=30, x=25, 
        null.value=0.6, alternative="greater", 
        succ.crit = "trial", alpha.final = 0.016,  
        a = 1, b=1) 

succ_ia_betabinom_one( N=40, n=30, x=25, 
        null.value=0.6, alternative="greater", test="exact",
        succ.crit = "trial", alpha.final = 0.016,  
        a = 1, b=1) 

succ_ia_betabinom_one( N=40, n=30, x=15, 
        null.value=0.6, alternative="greater", 
        succ.crit = "clinical", clin.succ.threshold =0.5,  
        a = 1, b=1) 
}

