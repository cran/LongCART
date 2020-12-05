\name{SurvCART}
\alias{SurvCART}
\title{Survival CART with time to event response via binary partitioning}
\description{Recursive partitioning for linear mixed effects model with survival 
data per SurvCART algorithm based on baseline partitioning variables (Kundu, 2020).
}
\usage{
SurvCART(data, patid, timevar, censorvar, gvars, tgvars, 
         time.dist="exponential", cens.dist="NA", event.ind=1, 
         alpha=0.05, minsplit=40, minbucket=20, print=FALSE)
}
\arguments{
\item{data}{name of the dataset. It must contain variable specified for \code{patid} (indicating subject id), all the variables specified in the formula and the baseline partitioning variables.}
\item{patid}{name of the subject id variable.}
\item{timevar}{name of the variable with follow-up times.}
\item{censorvar}{name of the variable with censoring status.}
\item{gvars}{list of partitioning variables of interest. Value of these variables should not change over time. Regarding categorical variables, only numerically coded categorical variables should be specified. For nominal categorical variables or factors, please first create corresponding dummy variable(s) and then pass through \code{gvars}. }
\item{tgvars}{types (categorical or continuous) of partitioning variables specified in \code{gvar}. For each of continuous partitioning variables, specify \code{1} and for each of the categorical partitioning variables, specify \code{0}. Length of \code{tgvars} should match to the length of \code{gvars}}
\item{time.dist}{name of time-to-event distribution. It can be one of the following distributions: \code{"exponential"}, \code{"weibull"}, \code{"lognormal"} or \code{"normal"}.}
\item{cens.dist}{name of censoring distribution. It can be one of the following distributions: \code{"exponential"}, \code{"weibull"}, \code{"lognormal"}, \code{"normal"} or \code{"NA"}. If specified \code{"NA"}, then parameter instability test corresponding to censoring distribution will not be performed.}
\item{event.ind}{value of the censoring variable indicating event.}
\item{alpha}{alpha (i.e., nominal type I error) level for parameter instability test}
\item{minsplit}{the minimum number of observations that must exist in a node in order for a split to be attempted.}
\item{minbucket}{he minimum number of observations in any terminal node.}
\item{print}{if \code{TRUE}, then summary such as number of subjects at risk, number of events, median event time and median censoring time model will be printed for each node.}
}
\details{
Construct survival tree based on heterogeneity in time-to-event and censoring distributions.


Exponential distribution: f(t)=lambda*exp(-lambda*t)

Weibull distribution: f(t)=alpha*lambda*t^(alpha-1)*exp(-lambda*t^alpha)

Lognormal distribution: f(t)=(1/t)*(1/sqrt(2*pi*sigma^2))*exp[-(1/2)*{(log(t)-mu)/sigma}^2]

Normal distribution: f(t)=(1/sqrt(2*pi*sigma^2))*exp[-(1/2)*{(t-mu)/sigma}^2]
}
\value{
\item{Treeout}{contains summary information of tree fitting for each terminal nodes and non-terminal nodes. Columns of \code{Treeout} include "ID", the (unique) node numbers that follow a binary ordering indexed by node depth, \code{n}, the number of subjectsreaching the node, \code{D}, the number of events reaching the node, \code{median.T}, the median survival time at the node, \code{median.C}, the median censoring time at the node, var, splitting variable, \code{index}, the cut-off value of splitting variable for binary partitioning, \code{p (Instability)}, the p-value for parameter instability test for the splitting variable, \code{loglik}, the log-likelihood of the node, \code{improve},  the improvement in deviance given by this split, and \code{Terminal},  indicator (True or False) of terminal node.  }
\item{logLik.tree}{log-likelihood of the tree-structured model, based on Cox model including sub-groups as covariates}
\item{logLik.root}{log-likelihood at the root node (i.e., without tree structure), based on Cox model without any covariate  }
\item{AIC.tree}{AIC of the tree-structured model, based on Cox model including sub-groups as covariates}
\item{AIC.root}{AIC at the root node (i.e., without tree structure), based on Cox model without any covariate }
\item{subj.class}{Assigned node for each individual subjects per fitted longitudinal tree}
\item{frame}{rpart compatible object}
\item{splits}{rpart compatible object}
\item{cptable}{rpart compatible object}
\item{functions}{rpart compatible object}
}
\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G.(2020). Survival trees based on score based parameter instability test. 
}
\seealso{
\code{\link{plot.SurvCART}}, \code{\link{KMPlot.SurvCART}}, \code{\link{text.SurvCART}}, \code{\link{StabCat.surv}}, \code{\link{StabCont.surv}}
}


\examples{
#--- Get the data
data(GBSG2)

#numeric coding of character variables
GBSG2$horTh1<- as.numeric(GBSG2$horTh)
GBSG2$tgrade1<- as.numeric(GBSG2$tgrade)
GBSG2$menostat1<- as.numeric(GBSG2$menostat)

#Add subject id
GBSG2$subjid<- 1:nrow(GBSG2)

#--- Run SurvCART() with time-to-event distribution: exponential, censoring distribution: None  
out<- SurvCART(data=GBSG2, patid="subjid", censorvar="cens", timevar="time", 
        gvars=c('horTh1', 'age', 'menostat1', 'tsize', 'tgrade1', 'pnodes', 'progrec', 'estrec'),  
        tgvars=c(0,1,0,1,0,1, 1,1),          
        event.ind=1,  alpha=0.05, minsplit=80, minbucket=40, print=TRUE)

#--- Plot tree
par(xpd = TRUE)
plot(out, compress = TRUE)
text(out, use.n = TRUE)

\donttest{
#Plot KM plot for sub-groups identified by tree
KMPlot.SurvCART(out, scale.time=365.25, type=1)


#--- Run SurvCART() with time-to-event distribution: weibull censoring distribution: None  
out2<- SurvCART(data=GBSG2, patid="subjid", censorvar="cens", timevar="time",  
        gvars=c('horTh1', 'age', 'menostat1', 'tsize', 'tgrade1', 'pnodes', 'progrec', 'estrec'),  
        tgvars=c(0,1,0,1,0,1, 1,1),          
        time.dist="weibull", event.ind=1, alpha=0.05, minsplit=80, minbucket=40, print=TRUE)


#--- Run SurvCART() with time-to-event distribution: weibull censoring distribution: exponential
out<- SurvCART(data=GBSG2, patid="subjid", censorvar="cens", timevar="time",  
        gvars=c('horTh1', 'age', 'menostat1', 'tsize', 'tgrade1', 'pnodes', 'progrec', 'estrec'),  
        tgvars=c(0,1,0,1,0,1, 1,1),          
        time.dist="weibull", cens.dist="exponential", event.ind=1, 
        alpha=0.05, minsplit=80, minbucket=40, print=TRUE)
}
}
