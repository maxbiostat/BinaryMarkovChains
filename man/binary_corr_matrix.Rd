% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary_multiESS.R
\name{binary_corr_matrix}
\alias{binary_corr_matrix}
\title{Correlation matrix for binary variables.}
\usage{
binary_corr_matrix(x, ncores = 2)
}
\arguments{
\item{x}{a matrix or data.frame containing only binary variables.}

\item{ncores}{(optional) number of cores to be used.}
}
\value{
a correlation matrix.
}
\description{
Correlation matrix for binary variables.
}
\examples{
N <- 100
M <- 1000
X <- matrix(rbinom(n = N*M, size = 1, p = .234), ncol = N, nrow = M)
cormat.binary <- binary_corr_matrix(X)
cormat.default <- cor(X)
plot(cormat.binary ~ cormat.default)
abline(a = 0, b = 1, lty = 2)
}
