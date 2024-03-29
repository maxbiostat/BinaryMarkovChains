% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stationarity.R
\name{compute_Hsquared}
\alias{compute_Hsquared}
\title{Compute the likelihood ratio statistic comparing a stationary Markov chain model for \code{seq} to a nonstationary model}
\usage{
compute_Hsquared(seq, L = 10, states = NULL, eps = 1e-17)
}
\arguments{
\item{seq}{a vector containing observations from an stochastic process.}

\item{L}{an integer specifying how many windows will be compared}

\item{states}{the possible states of the Markov chain. If \code{NULL},
will be deduced from 'seq'.}

\item{eps}{a small number to avoid zero probabilities}
}
\value{
a likelihood ratio statistic comparing H1 = not stationary vs
H0 = stationary.
}
\description{
Compute the likelihood ratio statistic comparing a stationary Markov chain model for \code{seq} to a nonstationary model
}
\examples{
X <- rbinom(n = 1E5, size = 1, prob = 0.8)
compute_Hsquared(X)
}
