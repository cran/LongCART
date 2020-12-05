\name{plot}
\alias{plot}
\alias{plot.SurvCART}
\alias{plot.LongCART}
\title{Plot an SurvCART or LongCART Object}
\description{
Plots an SurvCART or LongCART object on the current graphics device.
}
\usage{
\method{plot}{SurvCART}(x, uniform = FALSE, branch = 1, compress = FALSE, 
              nspace = branch, margin = 0, minbranch = 0.3, ...)
\method{plot}{LongCART}(x, uniform = FALSE, branch = 1, compress = FALSE, 
              nspace = branch, margin = 0, minbranch = 0.3, ...)
}
\arguments{
\item{x}{a fitted object of class \code{"SurvCART"}, containing a survival tree or \code{"LongCART"}, containing a longitudinal tree.}
\item{uniform}{similar to \code{plot.rpart}; if TRUE, uniform vertical spacing
 of the nodes is used; this may be less cluttered when fitting a large plot
onto a page. The default is to use a non-uniform spacing proportional to the
error in the fit.}
\item{branch}{similar to \code{plot.rpart};
controls the shape of the branches from parent to child node. Any number from
0 to 1 is allowed. A value of 1 gives square shouldered branches, a value of
0 give V shaped branches, with other values being intermediate.}
\item{compress}{similar to \code{plot.rpart};
if \code{FALSE}, the leaf nodes will be at the horizontal plot coordinates of \code{1:nleaves}. If \code{TRUE}, the routine attempts a more compact arrangement of the tree.}
\item{nspace}{similar to \code{plot.rpart}; the amount of extra space between a node with children
and a leaf, as compared to the minimal space between leaves. Applies to compressed trees only.
The default is the value of \code{branch}.}
\item{margin}{similar to \code{plot.rpart}; an extra fraction of white space to
leave around the borders of the tree. (Long labels sometimes get cut off by the default computation). }
\item{minbranch}{similar to \code{plot.rpart}; set the minimum length for a branch to \code{minbranch}
times the average branch length. This parameter is ignored if \code{uniform=TRUE}. Sometimes a split will
give very little improvement, or even (in the classification case) no improvement at all. A tree with
branch lengths strictly proportional to improvement leaves no room to squeeze in node labels.}
\item{...}{arguments to be passed to or from other methods.}
}
\details{
This function is a method for the generic function plot, for objects of class
\code{SurvCART}. The y-coordinate of the top node of the tree will always be 1.}
\value{
The coordinates of the nodes are returned as a list, with components \code{x}
 and \code{y}.
}
\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G.(2020). Survival trees based on score based parameter instability test. 
}
\seealso{
\code{\link{text}}, \code{\link{SurvCART}}, \code{\link{LongCART}}
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