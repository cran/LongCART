\name{ProfilePlot}
\alias{ProfilePlot}
\title{
Population level longitudinal profile plot for sub-groups
}
\description{
Generates population level longitudinal profile plot for each of sub-groups (i.e., terminal nodes) associated with longitudinal tree generated by LongCART()
}
\usage{
ProfilePlot(x, timevar, timevar.power=NULL, covariate.val=NULL,
                xlab=NULL, ylab=NULL, sg.title=4, mfrow=NULL, ...)
}
\arguments{
  \item{x}{a fitted object of class \code{"LongCART"}, containing a longitudinal tree.}
  \item{timevar}{Speciy the variable name contining time informaiton in the dataset that was used to fit LongCART object}
  \item{timevar.power}{Mandatory when the fixed part of the fitted model contains term as time with power not equal to  1. For example, if fixed part of the model is t + sqrt{t} + cov1, then specify \code{c(1, 0.5, NA}. If fixed part of the model is t + t^2 + cov1, then specify \code{c(1, 2, NA)}.}
  \item{covariate.val}{Specify the covariate values for generation of longitudinal profile plot. Iin the longitudinal profile plot, only time can vary and therefore, and therefore the value for the other covariates are fixed at constant value. This is not needed if the longitudinal model does not contain additional covariate(s). By default, the covariates values are specified at median value over all the datapoint (not at the subject level).  For example, if the fixed part of the model is t + cov1, then \code{c(NA, 100)} sets the value of cov1 at 100. Similarly, If fixed part of the model is t + t^2 + cov1, then \code{c(NA, NA, 100)} would be acceptable.}
  \item{xlab}{Optional label for X-axis}
  \item{ylab}{Optional label for Y-axis}
  \item{sg.title}{1 for sub-groups' title as \code{Sub-group=x}, 2 for sub-groups' title as \code{Node=x}, 3 for sub-groups' title as \code{Sub-group=x (Node=x)}, and 4 for sub-groups' title as \code{<node number>: <node defintion>}.}
  \item{mfrow}{Desired frame for fitting multiple plots. Default option is to include plots for all subgroups in the same frame. This  input is ignored when \code{overlay=TRUE}.}
  \item{...}{Graphical parameters other than \code{x, y, type, xlab, ylab}.}
}
\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G., and Harezlak, J. (2019). Regression trees for longitudinal data with baseline covariates. Biostatistics & Epidemiology, 3(1):1-22.
}
\seealso{
\code{\link{text}}, \code{\link{plot}}, \code{\link{LongCART}}
}

\examples{
\donttest{
#--- Get the data
data(ACTG175)

#-----------------------------------------------#
#   model: cd4~ time + subject(random)          #
#-----------------------------------------------#

#--- Run LongCART()
gvars=c("gender", "wtkg", "hemo", "homo", "drugs",
        "karnof", "oprior", "z30", "zprior", "race",
        "str2", "symptom", "treat", "offtrt")
tgvars=c(0, 1, 0, 0, 0,
         1, 0, 0, 0, 0,
         0, 0, 0, 0)


out1<- LongCART(data=ACTG175, patid="pidnum", fixed=cd4~time,
                gvars=gvars, tgvars=tgvars, alpha=0.05,
                minsplit=100, minbucket=50, coef.digits=2)

#--- Plot longitudinal profiles of subgroups
ProfilePlot(x=out1, timevar="time")

#-----------------------------------------------#
#   model: cd4~ time+ time^2 + subject(random)  #
#-----------------------------------------------#

ACTG175$time2<- ACTG175$time^2

out2<- LongCART(data=ACTG175, patid="pidnum", fixed=cd4~time + time2,
                gvars=gvars, tgvars=tgvars, alpha=0.05,
                minsplit=100, minbucket=50, coef.digits=2)

ProfilePlot(x=out2, timevar="time", timevar.power=c(1,2))


#--------------------------------------------------------#
#   model: cd4~ time+ time^2 + subject(random) + karnof  #
#--------------------------------------------------------#

out3<- LongCART(data=ACTG175, patid="pidnum", fixed=cd4~time + time2 + karnof,
                gvars=gvars, tgvars=tgvars, alpha=0.05,
                minsplit=100, minbucket=50, coef.digits=2)

#the value of the covariate karnof is set at median by default
ProfilePlot(x=out3, timevar="time", timevar.power=c(1,2, NA))

#the value of the covariate karnof is set at 120
ProfilePlot(x=out3, timevar="time", timevar.power=c(1,2, NA), 
                    covariate.val=c(NA, NA, 120))

}
}
