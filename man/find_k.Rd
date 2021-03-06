% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_selection.r
\name{find_k}
\alias{find_k}
\title{Thin a sample until the first-order Markov model is preferred.}
\usage{
find_k(seq, max_k = 10)
}
\arguments{
\item{seq}{a vector containing observations from an stochastic process.}

\item{max_k}{maximum thinning to be tried.}
}
\value{
a list containing the selected \code{k}, the achieved Bayes factor (\code{BF})
and whether one could indeed find a \code{k} less than \code{max_k}
such that the first-order model is preferred.
}
\description{
Thin a sample until the first-order Markov model is preferred.
}
\details{
This is the same approach as Raftery & Lewis (1992) "How many iterations in the Gibbs sampler?"
}
\examples{
library(markovchain)
M <- 1E4
p <-  runif(1)
( maxA <- min(p/(1-p), 1) ) ## maximum alpha
a <- runif(1, 0, maxA) 
b <- exp(log(a)  + log(1-p) - log(p)) 
mat <- matrix(c(1-a, a, b, 1-b), ncol = 2, nrow = 2, byrow = TRUE)
MC.binary <- new("markovchain",
 states = c("0", "1"),
  transitionMatrix = mat, name = "Binary")
out1 <- markovchainSequence(n = M, markovchain = MC.binary, t0 = "0")
X1 <- as.numeric(out1)

rs <- runif(4^2)
mat2 <- structure(c(0.264985809735216, 0.598430024309721, 0.140708602965477, 
0.0566167324620894, 0.278728488803952, 0.0427599923571986, 0.453529988769559, 
0.604004074666218, 0.116827824413052, 0.0358869202453851, 0.181735130309576, 
0.266807259734723, 0.33945787704778, 0.322923063087696, 0.224026277955388, 
0.0725719331369693), .Dim = c(4L, 4L))
 MC.binary2nd <- new("markovchain",
  states = c("00", "01", "10", "11"),
   transitionMatrix = mat2, name = "Binary2nd")
out2 <- markovchainSequence(n = M/2, markovchain = MC.binary2nd, t0 = "01")
X2 <- as.numeric(strsplit(paste(out2, collapse = ""), "")[[1]])
find_k(X1)
find_k(X2)
}
