\name{PoS}
\alias{PoS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability of trial ans clinical success for a prospective trial using normal-normal approximation
}
\description{
This function can be used to determine probability of trial success and clinical success based on the prior distribution for each of continuous, binary  and time-to-event endpoints. The calculation is carried out assuming normal distribution for estimated parameter and normal prior distribution.
}
\usage{
PoS(type, nsamples, null.value = NULL, alternative = "greater", 
    N = NULL, D = NULL, a = 1, 
    succ.crit = "trial", Z.crit.final = 1.96, alpha.final = 0.025, 
    clin.succ.threshold = NULL, se.exp = NULL, 
    meandiff.prior = NULL, mean.prior = NULL, sd.prior = NULL, 
    propdiff.prior = NULL, prop.prior = NULL, hr.prior = NULL, D.prior = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{Type of the endpoint. It could be \code{cont} for continuous, \code{bin} for binary and \code{surv} for survival endpoint.}
  \item{nsamples}{Number of samples. For continuous and binary case, it can be 1 or 2. For survival endpoint, it can be only 2. }
  \item{null.value}{The specified value under null hypothesis. Default is 0 for continuous and binomial case and 1 for survival case.}
  \item{alternative}{Direction of alternate hypothesis. Can be "greater" or "less". Default is "less" for test of HR and "greater" otherwise.}
  \item{N}{Total sample size at final analysis. Cannot be missing for continuous and binary endpoint.}
  \item{D}{Total number of events at final analysis. Cannot be missing for survival endpoint.}
  \item{a}{Allocation ratio in two sample case.}
  \item{succ.crit}{Specify "trial" for trial success (i.e., null hypothesis is rejected at final analysis) or "clinical" for clinical success (i.e., estimated value at the final analysis is greater than clinically meaningful value as specified under \code{clin.succ.threshold}.)}
  \item{Z.crit.final}{The rejection boundary at final analysis in Z-value scale. Either \code{alpha.final} or \code{Z.crit.final} must be specified when determining trial success.}
  \item{alpha.final}{The rejection boundary at final analysis in alpha (1-sided) scale (e.g., 0.025). Either \code{alpha.final} or \code{Z.crit.final} must be specified when determining trial success.}
  \item{clin.succ.threshold}{Clinically meaningful value. Required when \code{succ.crit="clinical"}.}
  \item{se.exp}{Expected standard error to be observed in the study. Must be specified in continuous case and two-sample binary case. }
  \item{meandiff.prior}{Mean value of prior distribution for mean difference. Relevant for two-sample continuous case.}
  \item{mean.prior}{Mean value of prior distribution for mean. Relevant for one-sample continuous case.}
  \item{sd.prior}{Standard deviation of prior distribution for mean difference (2-sample continuous case) or mean (1-sample continuous case) or prop (2-sample binary case) or difference of proportion (1-sample binary case) or log(HR) (2 sample survival case). }
  \item{propdiff.prior}{Mean value of prior distribution for difference in proportion. Relevant for two-sample binomial case.}
  \item{prop.prior}{Mean value of prior distribution for proportion. Relevant for one-sample binomial case.}
  \item{hr.prior}{Mean value of prior distribution for hazards ratio (HR). Relevant for two-sample survival case.}
  \item{D.prior}{Ignored if \code{sd.prior} is specified. If \code{sd.prior} is not specified then sd.prior is determined as 2/\code{D.prior}. Relevant for two-sample survival case.}
}

\details{
This function can be used to determine probability of success (PoS) for a prospective trial for each of continuous (one-sample or two-samples), binary (one-sample or two-samples) and time-to-event endpoints (two-samples). The PoS is calculated based on the prior distribution and expected standard error in the estimate in trial. The calculation PoS is carried out assuming normal distribution for estimated parameter and normal prior distribution. This function can be used to determine clinical success (\code{succ.crit="clinical"}) and trial success (\code{succ.crit="trial"}). For clinical success, \code{clin.succ.threshold} must be specified. For trial success, \code{Z.crit.final or alpha.final} must be specified.
}

\author{Madan Gopal Kundu <madan_g.kundu@yahoo.com>}

\seealso{\code{\link{succ_ia_betabinom_one}}, \code{\link{succ_ia_betabinom_two}}, \code{\link{succ_ia}}   }

\examples{
#--- Example 1
PoS(type="cont", nsamples=2, null.value=-0.05, alternative="greater", 
        N=1552, a=1,  
        succ.crit="trial", Z.crit.final=1.97,
        se.exp=0.12*sqrt(1/776 + 1/776),
        meandiff.prior=0, sd.prior=0.02) 

#--- Example 2
PoS(type="bin", nsamples=2, null.value=0, alternative="greater", 
        N=210, a=2,  
        succ.crit="trial", Z.crit.final=2.012,
        se.exp=0.5*sqrt(1/140 + 1/70),
        propdiff.prior=0.20, sd.prior=sqrt(0.06)) 

PoS(type="bin", nsamples=2, null.value=0, alternative="greater", 
        N=210, a=2,  
        succ.crit="clinical", clin.succ.threshold =0.15,
        se.exp=0.5*sqrt(1/140 + 1/70),
        propdiff.prior=0.20, sd.prior=sqrt(0.06)) 

#--- Example 4
PoS(type="surv", nsamples=2, null.value=1, alternative="less", 
        D=441,  
        succ.crit="trial", Z.crit.final=1.96,
        hr.prior=0.71, D.prior=133) 

PoS(type="surv", nsamples=2, null.value=1, alternative="less", 
        D=441,  
        succ.crit="clinical", clin.succ.threshold =0.8,
        hr.prior=0.71, D.prior=133) 

}
