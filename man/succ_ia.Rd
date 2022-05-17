\name{succ_ia}
\alias{succ_ia}
\title{
Conditional power and predictive power of success based on interim results using normal-normal approximation
}
\description{
This function can be used to determine conditional power and predictive power for trial success and clinical success based on the interim results and prior distribution for each of continuous, binary  and time-to-event endpoints. The calculation is carried out assuming normal distribution for estimated parameter and normal prior distribution.
}
\usage{
succ_ia(type, nsamples, null.value = NULL, alternative = NULL, 
        N = NULL, n = NULL, D = NULL, d = NULL, a = 1, 
        meandiff.ia = NULL, mean.ia = NULL, 
        propdiff.ia = NULL, prop.ia = NULL, hr.ia = NULL, 
        stderr.ia = NULL, sd.ia = NULL, 
        succ.crit = "trial", Z.crit.final = 1.96, 
        alpha.final = 0.025, clin.succ.threshold = NULL, 
        meandiff.exp = NULL, mean.exp = NULL, 
        propdiff.exp = NULL, prop.exp = NULL, hr.exp = NULL, 
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
  \item{n}{Total sample size at interim analysis. Cannot be missing for continuous and binary endpoint.}
  \item{D}{Total number of events at final analysis. Cannot be missing for survival endpoint.}
  \item{d}{Total number of events at interim analysis. Cannot be missing for survival endpoint.}
  \item{a}{Allocation ratio in two sample case.}
  \item{meandiff.ia}{Estimated mean difference at interim analysis. Mandatory for continuous two sample case.}
  \item{mean.ia}{Estimated mean value at interim analysis. Mandatory for continuous single sample case}
  \item{propdiff.ia}{Estimated difference in proportion at interim analysis. Mandatory for binary two sample case}
  \item{prop.ia}{Estimated proportion at interim analysis. Mandatory for binary single sample case}
  \item{hr.ia}{Estimate hazards ratio (HR) at interim analysis. Mandatory for continuous single sample case}
  \item{stderr.ia}{Standard error (SE) of estimated mean difference (in one-sample continuous  case) or estimated mean (in two-sample continuous case) or estimated difference in proportion (in two-sample binary case) at interim analysis. For continuous case, if not specified, then the function attempts to estimate SE  from \code{sd.ia}. Mandatory for two-sample binary case.}
  \item{sd.ia}{Standard deviation of estimated mean difference (in one-sample continuous  case) or estimated mean (in two-sample continuous case) at interim analysis. If \code{stderr.ia} is specified, then the value of \code{sd.ia} is ignored. If code{stderr.ia} is not specified, then mandatory for continuous case.}
  \item{succ.crit}{Specify "trial" for trial success (i.e., null hypothesis is rejected at final analysis) or "clinical" for clinical success (i.e., estimated value at the final analysis is greater than clinically meaningful value as specified under \code{clin.succ.threshold}.)}
  \item{Z.crit.final}{The rejection boundary at final analysis in Z-value scale. Either \code{alpha.final} or \code{Z.crit.final} must be specified when determining trial success.}
  \item{alpha.final}{The rejection boundary at final analysis in alpha (1-sided) scale (e.g., 0.025). Either \code{alpha.final} or \code{Z.crit.final} must be specified when determining trial success.}
  \item{clin.succ.threshold}{Clinically meaningful value. Required when \code{succ.crit="clinical"}.}
  \item{meandiff.exp}{Expected mean difference in post interim data. Relevant for two-sample continuous case.}
  \item{mean.exp}{Expected mean in post interim data. Relevant for one-sample continuous case.}
  \item{propdiff.exp}{Expected difference in proportion in post interim data. Relevant for two-sample binary case.}
  \item{prop.exp}{Expected proportion in post interim data. Relevant for one-sample binary case.}
  \item{hr.exp}{Expected hazards ratio (HR) in post interim data. Relevant for two-sample survival case.}
  \item{meandiff.prior}{Mean value of prior distribution for mean difference. Relevant for two-sample continuous case.}
  \item{mean.prior}{Mean value of prior distribution for mean. Relevant for one-sample continuous case.}
  \item{sd.prior}{Standard deviation of prior distribution for mean difference (2-sample continuous case) or mean (1-sample continuous case) or prop (2-sample binary case) or difference of proportion (1-sample binary case) or log(HR) (2 sample survival case). }
  \item{propdiff.prior}{Mean value of prior distribution for difference in proportion. Relevant for two-sample binomial case.}
  \item{prop.prior}{Mean value of prior distribution for proportion. Relevant for one-sample binomial case.}
  \item{hr.prior}{Mean value of prior distribution for hazards ratio (HR). Relevant for two-sample survival case.}
  \item{D.prior}{Ignored if \code{sd.prior} is specified. If \code{sd.prior} is not specified then sd.prior is determined as 2/\code{D.prior}. Relevant for two-sample survival case.}
}
\details{
This function can be used to determine Conditional power (CP) and Predictive power or predictive probability of success (PPoS) based on the interim results for each of continuous (one-sample or two-samples), binary (one-sample or two-samples) and time-to-event endpoints (two-samples). The PPoS can be  based on interim results only or using both prior information and interim results. The calculation of CP and PPoS are carried out assuming normal distribution for estimated parameter and normal prior distribution. This function can be used to determine clinical success (\code{succ.crit="clinical"}) and trial success (\code{succ.crit="trial"}). For clinical success, \code{clin.succ.threshold} must be specified. For trial success, \code{Z.crit.final or alpha.final} must be specified.

In order to calculate CP and PPoS, \code{succ.ia()} should be invoked in the following form:

Continuous-two sample case (trial success): 

\code{succ.ia(type="cont", nsamples=2, null.value=, alternative=,  N=, n=, a, meandiff.ia, stderr.ia=, succ.crit="trial", Z.crit.final=)}

Continuous-two sample case (clinical success): 

\code{succ.ia(type="cont", nsamples=2, null.value=, alternative=,  N=, n=, a, meandiff.ia, stderr.ia=, succ.crit="clinical", clin.succ.threshold=)}

Continuous-one sample case (trial success): 

\code{succ.ia(type="cont", nsamples=1, null.value=, alternative=,  N=, n=, mean.ia, stderr.ia=, succ.crit="trial", Z.crit.final=)}

Continuous-one sample case (clinical success): 

\code{succ.ia(type="cont", nsamples=1, null.value=, alternative=,  N=, n=, mean.ia, stderr.ia=, succ.crit="clinical", clin.succ.threshold=)}

Binary-two sample case (trial success): 

\code{succ.ia(type="bin", nsamples=2, null.value=, alternative=,  N=, n=, a, propdiff.ia, stderr.ia=, succ.crit="trial", Z.crit.final=)}

Binary-two sample case (clinical success): 

\code{succ.ia(type="bin", nsamples=2, null.value=, alternative=,  N=, n=, a, propdiff.ia, stderr.ia=, succ.crit="clinical", clin.succ.threshold=)}

Binary-one sample case (trial success): 

\code{succ.ia(type="bin", nsamples=1, null.value=, alternative=,  N=, n=, prop.ia, succ.crit="trial", Z.crit.final=)}

Binary-one sample case (clinical success): 

\code{succ.ia(type="bin", nsamples=1, null.value=, alternative=,  N=, n=, prop.ia, succ.crit="clinical", clin.succ.threshold=)}

Survival-two sample case (trial success): 

\code{succ.ia(type="surv", nsamples=2, null.value=, alternative=,  D=, d=, a, hr.ia, succ.crit="trial", Z.crit.final=)}

Survival-two sample case (clinical success): 

\code{succ.ia(type="surv", nsamples=2, null.value=, alternative=,  D=, d=, a, hr.ia, succ.crit="clinical", clin.succ.threshold=)}

The conditional power is calculated assuming interim trend for post-interim data. If \code{meandiff.exp} (for continuous 2-samples case), \code{mean.exp} (for continuous 1-sample case), \code{propdiff.exp} (for binomial 2-samples case), \code{prop.exp} (for binomial 1-sample case), or \code{hr.exp} (for survival 2-samples case) is specified, then conditional power would be calculated using these specified value as well.

The Predictive power or Predictive probability of success (PPoS) is calculated based interim results. On top of this, it can also incorporate prior information. The prior information can be specified as follows: If \code{meandiff.prior, sd.prior} for continuous 2-samples case, \code{mean.prior, sd.prior} for continuous 1-sample case, \code{propdiff.prior, sd.prior} for binomial 2-samples case, \code{prop.prior, sd.prior} for binomial 1-sample case, and \code{hr.exp, sd.prior} (or, \code{hr.exp, D.prior}) for survival 2-samples case. 
}

\author{Madan Gopal Kundu <madan_g.kundu@yahoo.com>}
\references{
Kundu, M. G., Samanta, S., and Mondal, S. (2021). An introduction to the determination of the probability of a successful trial: Frequentist and Bayesian approaches. arXiv preprint arXiv:2102.13550.
}
\seealso{\code{\link{succ_ia_betabinom_one}}, \code{\link{succ_ia_betabinom_two}}, \code{\link{PoS}}   }

\examples{
#--- Lan et al. (2009), see #6. Example, outcome: Matching
succ_ia(type="cont", nsamples=1, null.value=0, alternative="greater", 
        N=225, n=45,   
        mean.ia=0, stderr.ia=1,       
        succ.crit="trial", Z.crit.final=1.96) 


#--- Dallow et al. (2011), see Figure 1. Example, outcome: Matching
succ_ia(type="cont", nsamples=1, null.value=0, alternative="greater", 
        N=100, n=50,   
        mean.ia=1.364, stderr.ia=1,      
        succ.crit="trial",  Z.crit.final=1.64) 

#--- Example 1 in the paper (Continuous endpoint)
succ_ia(type="cont", nsamples=2, null.value=-0.05, alternative="greater",
        N=1552, n=776, a=1,   
        meandiff.ia=-0.025, sd.ia=0.16,      
        succ.crit="trial", Z.crit.final=1.97,  
        meandiff.exp=-0.030,
        meandiff.prior=0, sd.prior=0.02) 

#--- Example 2 in the paper (Binary endpoint)
p1<- 0.379; p2<- 0.222
n1<- 105; n2<- 53

#-- Trial success
succ_ia(type="bin", nsamples=2, null.value=0, alternative="greater",
        N=210, n=158,  a=2,
        propdiff.ia=p1-p2,
        stderr.ia=sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2), 
        succ.crit="trial", Z.crit.final=2.012,
        propdiff.exp=0.20,
        propdiff.prior=0.20, sd.prior=sqrt(0.06))  

#-- Clinical success
succ_ia(type="bin", nsamples=2, null.value=0, alternative="greater",
        N=210, n=158,  a=2,
        propdiff.ia=p1-p2,
        stderr.ia=sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2), 
        succ.crit="clinical", clin.succ.threshold=0.15,
        propdiff.exp=0.20,
        propdiff.prior=0.20, sd.prior=sqrt(0.06)) 


#--- Example 3 in the paper (Survival endpoint)

#--- Trial success
succ_ia(type="surv", nsamples=2, null.value=1, alternative="less", 
        D=441, d=346, a=1,   
        hr.ia=0.82,        
        succ.crit="trial", Z.crit.final=2.012,            
        hr.exp=0.75,
        hr.prior=0.71, D.prior=133) 

#--- clinical success
succ_ia(type="surv", nsamples=2, null.value=1, alternative="less", 
        D=441, d=346, a=1,   
        hr.ia=0.82,        
        succ.crit="clinical", clin.succ.threshold=0.80,            
        hr.exp=0.75,
        hr.prior=0.71, D.prior=133) 

}
