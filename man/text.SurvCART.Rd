\name{text.SurvCART}
\alias{text.SurvCART}
\alias{text}
\title{Place text on SurvCART tree}
\description{
Labels the current plot of the tree generated from SurvCART object with text.
}
\usage{
\method{text}{SurvCART}(x, splits = TRUE, all = FALSE,
             use.n = FALSE, minlength = 1L, ...)
}
\arguments{
\item{x}{a fitted object of class \code{"SurvCART"}, containing a survival tree.}
\item{splits}{similar to text.rpart; 
logical flag. If \code{TRUE} (default), then the splits in the tree are labeled with the criterion for the split.}
\item{all}{similar to text.rpart; 
Logical. If TRUE, all nodes are labeled, otherwise just terminal nodes.}
\item{use.n}{similar to text.rpart; 
Logical. If TRUE, adds \code{n} to label.}
\item{minlength}{similar to text.rpart; 
the length to use for factor labels. A value of 1 causes them to be printed as 'a', 'b', ..... Larger values use abbreviations of the label names. See the labels.rpart function for details.}
\item{...}{arguments to be passed to or from other methods.}
}

\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G.(2020). Survival trees based on score based parameter instability test. 
}
\seealso{
\code{\link{plot.SurvCART}}, \code{\link{SurvCART}}, \code{\link{KMPlot.SurvCART}}, \code{\link{StabCat.surv}}, \code{\link{StabCont.surv}}
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

#--- Run SurvCART()
out<- SurvCART(data=GBSG2, patid="subjid", censorvar="cens", timevar="time", event.ind=1, 
        gvars=c('horTh1', 'age', 'menostat1', 'tsize', 'tgrade1', 'pnodes', 'progrec', 'estrec'),  
        tgvars=c(0,1,0,1,0,1, 1,1),          
        alpha=0.05, minsplit=80,  
        minbucket=40, print=TRUE)

#--- Plot tree
par(xpd = TRUE)
plot(out, compress = TRUE)
text(out, use.n = TRUE)

}
