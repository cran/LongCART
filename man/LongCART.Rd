\name{LongCART}
\alias{LongCART}
\title{Longitudinal CART with continuous response via binary partitioning}
\description{Recursive partitioning for linear mixed effects model with continuous
univariate response variables per LonCART algorithm based on baseline partitioning
variables (Kundu and Harezlak, 2019).
}
\usage{
LongCART(data, patid, fixed, gvars, tgvars, minsplit=40,
         minbucket=20, alpha=0.05, coef.digits=2, print.lme=FALSE)
}
\arguments{
\item{data}{name of the dataset. It must contain variable specified for \code{patid} (indicating subject id), all the variables specified in the formula and the baseline partitioning variables.}
\item{patid}{name of the subject id variable.}
\item{fixed}{a two-sided linear formula object describing the fixed-effects part of the model, with the response on the left of a \code{~ }operator and the terms, separated by \code{+} operators, on the right. Model with \code{-1} to the end of right side indicates no intercept. For model with no fixed effect beyond intercept, please specify only \code{1} right to the \code{~} operator.}
\item{gvars}{list of partitioning variables of interest. Value of these variables should not change over time. Regarding categorical variables, only numerically coded categorical variables should be specified. For nominal categorical variables or factors, please first create corresponding dummy variable(s) and then pass through \code{gvars}. }
\item{tgvars}{types (categorical or continuous) of partitioning variables specified in \code{gvar}. For each of continuous partitioning variables, specify \code{1} and for each of the categorical partitioning variables, specify \code{0}. Length of \code{tgvars} should match to the length of \code{gvars}}
\item{minsplit}{the minimum number of observations that must exist in a node in order for a split to be attempted.}
\item{minbucket}{he minimum number of observations in any terminal node.}
\item{alpha}{alpha (i.e., nominal type I error) level for parameter instability test}
\item{coef.digits}{decimal points for displaying coefficients in the tree  structure.}
\item{print.lme}{if \code{TRUE}, then summary of fitte model from \code{lme()} will be printed for each node.}
}
\details{
Construct regression tree based on heterogeneity in linear mixed effects models of following type:
\code{Y_i(t)= W_i(t)theta + b_i + epsilon_{it}}
where \code{W_i(t)} is the design matrix, \code{theta} is the parameter associated with \code{W_i(t)} and \code{b_i} is the random intercept.
Also, \code{epsilon_{it} ~  N(0,sigma ^2)} and \code{b_i ~ N(0, sigma_u^2)}.
}
\value{
\item{Treeout}{contains summary information of tree fitting for each terminal nodes and non-terminal nodes. Columns of \code{Treeout} include "ID", the (unique) node numbers that follow a binary ordering indexed by node depth, \code{n}, the number of observations reaching the node, \code{yval}, the fitted model of the response at the node, var, a factor giving the names of the variables used in the split at each, \code{index}, the cut-off value of splitting variable for binary partitioning, \code{p (Instability)}, the p-value for parameter instability test for the splitting variable, \code{loglik}, the log-likelihood of the node, \code{improve},  the improvement in deviance given by this split, and \code{Terminal},  indicator (True or False) of terminal node.  }
\item{p}{number of fixed parameters}
\item{AIC.tree}{AIC of the tree-structured model}
\item{AIC.root}{AIC at the root node (i.e., without tree structure)}
\item{improve.AIC}{improvement in AIC due to tree structure (AIC.tree - AIC.root)}
\item{logLik.tree}{log-likelihood of the tree-structured model}
\item{logLik.root}{log-likelihood at the root node (i.e., without tree structure) }
\item{Deviance}{2*(logLik.tree-logLik.root)}
\item{LRT.df}{degrees of freedom for likelihood ratio test comparing tree-structured model with the model at root node.}
\item{LRT.p}{p-value for likelihood ratio test comparing tree-structured model with the model at root node.}
\item{nodelab}{List of subgroups or terminal nodes with their description}
\item{varnam}{List of splitting variables}
\item{data}{the dataset originally supplied}
\item{patid}{the patid variable originally supplied}
\item{fixed}{the fixed part of the model originally supplied}
\item{frame}{rpart compatible object}
\item{splits}{rpart compatible object}
\item{cptable}{rpart compatible object}
\item{functions}{rpart compatible object}
}
\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G., and Harezlak, J. (2019). Regression trees for longitudinal data with baseline covariates. Biostatistics & Epidemiology, 3(1):1-22.
}
\seealso{
\code{\link{plot}}, \code{\link{text}}, \code{\link{ProfilePlot}}, \code{\link{StabCat}}, \code{\link{StabCont}}, \code{\link{predict}}
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

#--- Plot tree
par(mfrow=c(1,1))
par(xpd = TRUE)
plot(out1, compress = TRUE)
text(out1, use.n = TRUE)

#--- Plot longitudinal profiles of subgroups
ProfilePlot(x=out1, timevar="time")

#-----------------------------------------------#
#   model: cd4~ time+ time^2 + subject(random)  #
#-----------------------------------------------#

ACTG175$time2<- ACTG175$time^2

out2<- LongCART(data=ACTG175, patid="pidnum", fixed=cd4~time + time2,
                gvars=gvars, tgvars=tgvars, alpha=0.05,
                minsplit=100, minbucket=50, coef.digits=2)


par(mfrow=c(1,1))
par(xpd = TRUE)
plot(out2, compress = TRUE)
text(out2, use.n = TRUE)

ProfilePlot(x=out2, timevar="time", timevar.power=c(1,2))


#--------------------------------------------------------#
#   model: cd4~ time+ time^2 + subject(random) + karnof  #
#--------------------------------------------------------#

out3<- LongCART(data=ACTG175, patid="pidnum", fixed=cd4~time + time2 + karnof,
                gvars=gvars, tgvars=tgvars, alpha=0.05,
                minsplit=100, minbucket=50, coef.digits=2)


par(mfrow=c(1,1))
par(xpd = TRUE)
plot(out3, compress = TRUE)
text(out3, use.n = TRUE)

#the value of the covariate karnof is set at median by default
ProfilePlot(x=out3, timevar="time", timevar.power=c(1,2, NA)) 

#the value of the covariate karnof is set at 120
ProfilePlot(x=out3, timevar="time", timevar.power=c(1,2, NA), 
                     covariate.val=c(NA, NA, 120)) 

}
}
