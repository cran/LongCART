\name{GBSG2}
\alias{GBSG2}
\docType{data}
\title{
German Breast Cancer Study Group 2 (source: TH.data package)
}
\description{
A data frame containing the observations from the GBSG2 study.
}
\usage{data(GBSG2)}
\format{
  A data frame with 686 observations on the following 10 variables.
  \describe{
    \item{\code{horTh}}{hormonal therapy, a factor with levels \code{no} \code{yes}}
    \item{\code{age}}{age in years}
    \item{\code{menostat}}{menopausal status, a factor with levels \code{Pre} \code{Post}}
    \item{\code{tsize}}{tumor size (in mm)}
    \item{\code{tgrade}}{an ordered factor with levels \code{I} < \code{II} < \code{III}}
    \item{\code{pnodes}}{number of positive nodes}
    \item{\code{progrec}}{progesterone receptor (in fmol).}
    \item{\code{estrec}}{estrogen receptor (in fmol).}
    \item{\code{time}}{recurrence free survival time (in days).}
    \item{\code{cens}}{censoring indicator (0- censored, 1- event).}
  }
}

\references{
Schumacher M, Bastert G, Bojar H, Huebner K, Olschewski M, Sauerbrei W, Schmoor C, Beyerle C, Neumann RL, Rauschecker HF. Randomized 2 x 2 trial evaluating hormonal treatment and the duration of chemotherapy in node-positive breast cancer patients. German Breast Cancer Study Group. Journal of Clinical Oncology. 1994 Oct;12(10):2086-93.
}
\examples{
data(GBSG2)
}

