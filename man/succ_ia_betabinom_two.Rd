\name{succ_ia_betabinom_two}
\alias{succ_ia_betabinom_two}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determines predictive power of success based on interim results and beta priors for two-sample binary data
}
\description{
This function can be used to determine predictive power for trial success and clinical success based on the interim results and beta prior distribution for test of difference of two proportions. 
}
\usage{
succ_ia_betabinom_two(N.trt, N.con, 
                      n.trt, x.trt, n.con, x.con, 
                      alternative = "greater", test = "z", 
                      succ.crit = "trial", Z.crit.final = 1.96, 
                      alpha.final = 0.025, clin.succ.threshold = NULL, 
                      a.trt = 1, b.trt = 1, a.con = 1, b.con = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N.trt}{Sample size in treatment arm at final analysis. Cannot be missing.}
  \item{N.con}{Sample size in control arm at final analysis. Cannot be missing.}
  \item{n.trt}{Sample size in treatment arm at interim analysis. Cannot be missing.}
  \item{x.trt}{Number of observed response in treatment arm at interim analysis. Cannot be missing.}
  \item{n.con}{Sample size in control arm at interim analysis. Cannot be missing.}
  \item{x.con}{Number of observed response in control arm at interim analysis. Cannot be missing.}
  \item{alternative}{Direction of alternate hypothesis. Can be "greater" or "less".}
  \item{test}{Statistical test. Default is \code{"z"} for Z test. For Fisher's exact test, specify \code{"fisher"}.}
  \item{succ.crit}{Specify "trial" for trial success (i.e., null hypothesis is rejected at final analysis) or "clinical" for clinical success (i.e., estimated value at the final analysis is greater than clinically meaningful value as specified under \code{clin.succ.threshold}.)}
  \item{Z.crit.final}{The rejection boundary at final analysis in Z-value scale. Either \code{alpha.final} or \code{Z.crit.final} must be specified when determining trial success.}
  \item{alpha.final}{The rejection boundary at final analysis in alpha (1-sided) scale (e.g., 0.025). Either \code{alpha.final} or \code{Z.crit.final} must be specified when determining trial success.}
  \item{clin.succ.threshold}{Clinically meaningful value. Required when \code{succ.crit="clinical"}.}
  \item{a.trt}{Value of \code{a} corresponding to Beta(a, b) prior in treatment arm.}
  \item{b.trt}{Value of \code{b} corresponding to Beta(a, b) prior in treatment arm.}
  \item{a.con}{Value of \code{a} corresponding to Beta(a, b) prior for proportion in control arm.}
  \item{b.con}{Value of \code{b} corresponding to Beta(a, b) prior for proportion in control arm.}
}

\details{
This function can be used to determine Predictive power or predictive probability of success (PPoS) based on the interim results for comparison of two proportions. The calculation of PoS is carried out assuming beta prior distributions for proportions in both treatment and control arms. This function can be used to determine clinical success (\code{succ.crit="clinical"}) and trial success (\code{succ.crit="trial"}). For clinical success, \code{clin.succ.threshold} must be specified. For trial success, \code{Z.crit.final or alpha.final} must be specified.
}

\author{Madan Gopal Kundu <madan_g.kundu@yahoo.com>}

\seealso{\code{\link{succ_ia_betabinom_one}}, \code{\link{succ_ia}}, \code{\link{PoS}}    }

\examples{
\donttest{
succ_ia_betabinom_two( N.con=40, N.trt=40, 
        n.trt=30, x.trt=20, n.con=30, x.con=15, 
        alternative="greater", test="fisher",
        succ.crit = "trial", Z.crit.final = 1.96,  
        a.trt = 1, b.trt=1, a.con=1, b.con=1) 

succ_ia_betabinom_two( N.con=40, N.trt=40, 
        n.trt=30, x.trt=20, n.con=30, x.con=15, 
        alternative="greater", test="z",
        succ.crit = "trial", Z.crit.final = 1.96,  
        a.trt = 1, b.trt=1, a.con=1, b.con=1) 

succ_ia_betabinom_two( N.con=40, N.trt=40, 
        n.trt=30, x.trt=20, n.con=30, x.con=15, 
        alternative="greater", test="fisher",
        succ.crit = "clinical", clin.succ.threshold = 0.5, 
        a.trt = 1, b.trt=1, a.con=1, b.con=1) 

#--- Johns & Andersen, 1999, Example 1a (results matching)
succ_ia_betabinom_two( N.trt=32, N.con=32,  
        n.trt=12, x.trt=8, n.con=12, x.con=8, 
        alternative="greater", test="fisher",
        succ.crit = "clinical", clin.succ.threshold = 0,  
        a.trt = 1, b.trt=1, a.con=1, b.con=1) 

#--- Johns & Andersen, 1999, Example 1b (results matching)
succ_ia_betabinom_two( N.trt=32, N.con=32,  
        n.trt=12, x.trt=8, n.con=12, x.con=11, 
        alternative="greater", test="fisher",
        succ.crit = "clinical", clin.succ.threshold = 0,  
        a.trt = 1, b.trt=1, a.con=1, b.con=1) 

#--- Johns & Andersen, 1999, Example 2 (not matching, reported 0.586, got 0.536)
succ_ia_betabinom_two( N.trt=155+170, N.con=152+171,  
        n.trt=155, x.trt=13, n.con=152, x.con=21, 
        alternative="less", test="z",
        succ.crit = "trial", Z.crit.final = 1.96,  
        a.trt = 1, b.trt=1, a.con=1, b.con=1) 

succ_ia_betabinom_two( N.trt=155+170, N.con=152+171,  
        n.trt=155, x.trt=13, n.con=152, x.con=21, 
        alternative="less", test="fisher",
        succ.crit = "trial", Z.crit.final = 1.96,  
        a.trt = 1, b.trt=1, a.con=1, b.con=1) 
}
}
