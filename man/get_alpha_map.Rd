% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alpha_map.R
\name{get_alpha_map}
\alias{get_alpha_map}
\title{Maximum \emph{a posteriori} estimate of the alpha transition probability}
\usage{
get_alpha_map(dmat, k = 1, v = 1, p)
}
\arguments{
\item{dmat}{a 2 x 2 matrix with observed transitions.}

\item{k}{shape1 parameter of Beta prior on alpha (k > 0). Default is \code{k = 1}.}

\item{v}{shape2 parameter of Beta prior on alpha (v > 0). Default is \code{v = 1}.}

\item{p}{marginal success probability.}
}
\value{
MAP estimate of alpha
}
\description{
Maximum \emph{a posteriori} estimate of the alpha transition probability
}
\examples{
library(markovchain)
X <- sample(c(0, 1), 10000, replace = TRUE)
obs.trans <- createSequenceMatrix(X, sanitize = FALSE, possibleStates = c("0", "1"))
get_alpha_map(dmat = obs.trans, p = 1/2) ## should be close to 1/2 because of independent sampling
}
