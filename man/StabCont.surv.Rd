\name{StabCont.surv }
\alias{StabCont.surv }
\title{parameter stability test for continuous partitioning variable}
\description{Performs parameter stability test (Kundu, 2020) 
with continuous partitioning variable to determine whether the parameters 
of exponential time-to-event distribution and exponential censoring distribution remain 
same across all distinct values of given continupus partitioning variable.
}
\usage{
StabCont.surv(data, timevar, censorvar, splitvar, 
              time.dist="exponential", cens.dist="NA", event.ind=1, print=FALSE)
}
\arguments{
\item{data}{name of the dataset. It must contain variable specified for \code{timevar} (indicating follow-up times), \code{censorvar} (indicating censoring status) and the caterogrical partitioning variable of interest specified in \code{splitvar}. Note that, only numerically coded categorical variable should be specified. }
\item{timevar}{name of the variable with follow-up times.}
\item{censorvar}{name of the variable with censoring status.}
\item{time.dist}{name of time-to-event distribution. It can be one of the following distributions: \code{"exponential"}, \code{"weibull"}, \code{"lognormal"} or \code{"normal"}.}
\item{cens.dist}{name of censoring distribution. It can be one of the following distributions: \code{"exponential"}, \code{"weibull"}, \code{"lognormal"}, \code{"normal"} or \code{"NA"}. If specified \code{"NA"}, then parameter instability test corresponding to censoring distribution will not be performed.}
\item{event.ind}{value of the censoring variable indicating event.}
\item{splitvar}{the continuous partitioning variable of interest. }
\item{print}{if \code{TRUE}, then additional information including estimated parameters, score function and its variance will be printed.}
}
\details{
\code{StabCont.surv()} performs the following omnibus test

\code{H_0:theta_{(g)}=theta_0} vs. \code{H_1: theta_{(g)} ^= theta_0},     for all g

where, \code{theta_{(g)}} is the true value of \code{theta} for subjects with \code{X=C_g}. \code{theta} includes all the parameters of time to event distribution and also parameters of censoring distribution, if specified. \code{C_g} is the any value realized by continuous partitioning variable \code{X}.

Exponential distribution: f(t)=lambda*exp(-lambda*t)

Weibull distribution: f(t)=alpha*lambda*t^(alpha-1)*exp(-lambda*t^alpha)

Lognormal distribution: f(t)=(1/t)*(1/sqrt(2*pi*sigma^2))*exp[-(1/2)*{(log(t)-mu)/sigma}^2]

Normal distribution: f(t)=(1/sqrt(2*pi*sigma^2))*exp[-(1/2)*{(t-mu)/sigma}^2]
}
\value{
\item{pval}{p-value for parameter instability test}
\item{type}{1, if event times are more heterogeneous; 2, if censoring times are more hetergeneous.}
}
\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G., and Ghosh, S. (2021). Survival trees based on heterogeneity in time-to-event and censoring distributions using parameter instability test. Statistical Analysis and Data Mining: The ASA Data Science Journal, 14(5), 466-483.
}
\seealso{
\code{\link{StabCont.surv}}, \code{\link{SurvCART}}, \code{\link{plot}}, \code{\link{text}}
}

\examples{   
#--- time-to-event distribution: exponential, censoring distribution: None    
out1<- StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", event.ind=2) 
out1$pval

#--- time-to-event distribution: weibull, censoring distribution: None  
StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", 
              time.dist="weibull", event.ind=2) 

#--- time-to-event distribution: weibull, censoring distribution: exponential
StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", 
              time.dist="weibull", cens.dist="exponential", event.ind=2) 
}