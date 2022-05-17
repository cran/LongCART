\name{predict}
\alias{predict}
\alias{predict.SurvCART}
\alias{predict.LongCART}
\title{Predicts according to the fitted SurvCART or LongCART tree}
\description{
Predicts according to the fitted SurvCART or LongCART tree.
}
\usage{
\method{predict}{SurvCART}(object, newdata, ...)
\method{predict}{LongCART}(object, newdata, patid, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"SurvCART"}, containing a survival tree, or class \code{"LongCART"}, containing a longitudinal tree.}
\item{newdata}{The dataset for prediction.}
\item{patid}{Variable name containing patient id in the new dataset. Must for prediction based on LongCART object}
\item{...}{Please disregard.}
}
\details{
For prediction based on \code{"SurvCART"} algorithm, the predicted dataset includes the terminal node id the observation belongs to, and the median event and censoring times of the terminal id. 

For prediction based on \code{"LongCART"} algorithm, the predicted dataset includes the terminal node id the observation belongs to,  the fitted profile, and the predicted value based on the fitted profile. Note that the predicted value does not consider the random effects.  
}
\value{
For prediction based on \code{"SurvCART"} algorithm, the dataset adds to the following variables in the new dataset:
\item{node}{Terminal node id the observation belongs to}
\item{median.T}{Median event time of the terminal node id the observation belongs to}
\item{median.C}{Median censoring time of the terminal node id the observation belongs to}
\item{Q1.T}{First quartile for event time of the terminal node id the observation belongs to}
\item{Q1.C}{First quartile for censoring time of the terminal node id the observation belongs to}
\item{Q3.T}{Third quartile for event time of the terminal node id the observation belongs to}
\item{Q3.C}{Third quartile for censoring time of the terminal node id the observation belongs to}

For prediction based on \code{LongCART} algorithm, the dataset adds to the following variables in the new dataset:
\item{node.id}{Terminal node id the observation belongs to}
\item{profile}{The fitted profile of the terminal node id the observation belongs to}
\item{predval}{predicted value based on the fitted profile \code{profile}}
}
\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G., and Harezlak, J. (2019). Regression trees for longitudinal data with baseline covariates. Biostatistics & Epidemiology, 3(1):1-22.

Kundu, M. G., and Ghosh, S. (2021). Survival trees based on heterogeneity in time-to-event and censoring distributions using parameter instability test. Statistical Analysis and Data Mining: The ASA Data Science Journal, 14(5), 466-483.
}
\seealso{
\code{\link{SurvCART}}, \code{\link{LongCART}}
}

\examples{
\donttest{
#--- LongCART example

data(ACTG175)
gvars=c("gender", "wtkg", "hemo", "homo", "drugs",
        "karnof", "oprior", "z30", "zprior", "race",
        "str2", "symptom", "treat", "offtrt")
tgvars=c(0, 1, 0, 0, 0,
         1, 0, 0, 0, 0,
         0, 0, 0, 0)
out1<- LongCART(data=ACTG175, patid="pidnum", fixed=cd4~time,
                gvars=gvars, tgvars=tgvars, alpha=0.05,
                minsplit=100, minbucket=50, coef.digits=2)
pred1<- predict.LongCART(object=out1, newdata=ACTG175, patid="pidnum")
head(pred1)

#--- SurvCART example

data(GBSG2)
GBSG2$horTh1<- as.numeric(GBSG2$horTh)
GBSG2$tgrade1<- as.numeric(GBSG2$tgrade)
GBSG2$menostat1<- as.numeric(GBSG2$menostat)

GBSG2$subjid<- 1:nrow(GBSG2)

fit<- SurvCART(data=GBSG2, patid="subjid", censorvar="cens", timevar="time", 
        gvars=c('horTh1', 'age', 'menostat1', 'tsize', 'tgrade1', 'pnodes', 'progrec', 'estrec'),  
        tgvars=c(0,1,0,1,0,1, 1,1),          
        event.ind=1,  alpha=0.05, minsplit=80, minbucket=40, print=TRUE)

pred2<- predict.SurvCART(object=fit, newdata=GBSG2)
head(pred2)
}
}
