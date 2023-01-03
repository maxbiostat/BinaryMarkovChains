#' Compute transitions of a given order
#'
#' @param seq  a vector containing observations from an stochastic process.
#' @param order the order of the transitions to be computed (an integer).
#'
#' @return an array with the transitions in the format n_{ijk}, where i, j, k are in states(seq).
#' @export get_contingency_array
#'
#' @examples
#' X <- rbinom(n = 1E5, size = 1, prob = 0.8)
#' get_contingency_array(X, order = 3)
get_contingency_array <- function(seq, order, states = NULL){
  ## Many thanks to https://stackoverflow.com/users/3358272/r2evans
  ## https://stackoverflow.com/questions/67004206/speeding-up-r-code-to-compute-higher-order-transitions-in-a-markov-chain
  N <- length(seq)
  if(is.null(states)) states <- sort(unique(seq))
  nstates <- length(states)
  inds <- seq_along(states)
  K <- order + 1
  out <- array(0, dim = rep(nstates, K))
  for (z in K:N){
    pos <- matrix(match(seq[(z - K + 1):z], states), nrow = 1)
    out[pos] <- out[pos] + 1
  }
  return(out)
}
#' Compute the likelihood ratio statistic comparing a second-order to a first-order Markov chain model for \code{seq}
#'
#' @param seq a vector containing observations from an stochastic process.
#'
#' @return the likelihood ratio statistic (G^2) comparing a second-order model to a first-order one. 
#' @export compute_Gsquared
#' @details  Under the null, G^2 has a chi-square(2) distribution.
#' @examples
#' X <- rbinom(n = 1E5, size = 1, prob = 0.8)
#' compute_Gsquared(X)
compute_Gsquared <- function(seq){
  states <- sort(unique(seq))
  nstates <- length(states)
  dt <-  get_contingency_array(seq, 2)
  stats <- c()
  for(k in 1:nstates){
    for( i in 1:nstates){
      for(j in 1:nstates){
        nijk <-  dt[i , j, k] + 1E-17
        nijp <- sum(dt[i, j, ]) + 1E-17
        npjk <- sum(dt[, j, k]) + 1E-17
        npjp <- sum(dt[, j, ]) + 1E-17
        s <- nijk * log((nijk/nijp)/(npjk/npjp))
        stats <- c(stats, s)
      }
    }
  }
  out <- 2*sum(stats)
  return(out)
}
#' Compute the Bayes factor comparing a second-order Markov model to a first-order one.
#'
#' @param seq a vector containing observations from an stochastic process.
#'
#' @return the approximate Bayes factor comparing a second- to a first-order Markov chain model for \code{seq}.
#' @export compute_BF
#'
#' @examples
#' X <- rbinom(n = 1E5, size = 1, prob = 0.8)
#' compute_BF(X)
compute_BF <- function(seq){
  n <- length(seq)
  m <- length(unique(seq))
  Gsq <- compute_Gsquared(seq)
  BF <- (Gsq - 2*log(n))/2
  return(BF)
}
#' Thin a sample until the first-order Markov model is preferred.
#'
#' @param seq a vector containing observations from an stochastic process.
#' @param max_k maximum thinning to be tried.
#'
#' @return a list containing the selected \code{k}, the achieved Bayes factor (\code{BF})
#'  and whether one could indeed find a \code{k} less than \code{max_k}
#'   such that the first-order model is preferred.
#' @export find_k
#' @details This is the same approach as Raftery & Lewis (1992) "How many iterations in the Gibbs sampler?"
#' @examples
#' library(markovchain)
#' M <- 1E4
#' p <-  runif(1)
#' ( maxA <- min(p/(1-p), 1) ) ## maximum alpha
#' a <- runif(1, 0, maxA) 
#' b <- exp(log(a)  + log(1-p) - log(p)) 
#' mat <- matrix(c(1-a, a, b, 1-b), ncol = 2, nrow = 2, byrow = TRUE)
#' MC.binary <- new("markovchain",
#'  states = c("0", "1"),
#'   transitionMatrix = mat, name = "Binary")
#' out1 <- markovchainSequence(n = M, markovchain = MC.binary, t0 = "0")
#' X1 <- as.numeric(out1)
#' 
#' rs <- runif(4^2)
#' mat2 <- structure(c(0.264985809735216, 0.598430024309721, 0.140708602965477, 
#' 0.0566167324620894, 0.278728488803952, 0.0427599923571986, 0.453529988769559, 
#' 0.604004074666218, 0.116827824413052, 0.0358869202453851, 0.181735130309576, 
#' 0.266807259734723, 0.33945787704778, 0.322923063087696, 0.224026277955388, 
#' 0.0725719331369693), .Dim = c(4L, 4L))
#'  MC.binary2nd <- new("markovchain",
#'   states = c("00", "01", "10", "11"),
#'    transitionMatrix = mat2, name = "Binary2nd")
#' out2 <- markovchainSequence(n = M/2, markovchain = MC.binary2nd, t0 = "01")
#' X2 <- as.numeric(strsplit(paste(out2, collapse = ""), "")[[1]])
#' find_k(X1)
#' find_k(X2)
find_k <- function(seq, max_k = 10){
  k <- 1
  BF <- compute_BF(seq)
  while(BF > 0){
    k <- k + 1
    inds <- seq.int(1L, length(seq), k)
    BF <- compute_BF(seq = seq[inds])
    if(k >= max_k) break
  }
  return(
    list(
      k = k,
      BF = BF,
      success = (BF < 0) 
    )
  )
}