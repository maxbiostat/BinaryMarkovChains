#' Return the frequencies of each state in each window.
#'
#' @param samples a collection of N samples from a two-state Markov chain.
#' @param  k a number of windows into which to break up the samples.
#'
#' @return a matrix with counts of each state in each window.
#' @export split_counts
#'
#' @examples
#' X <- rbinom(1000, size = 1, prob =.8)
#' tab <- split_counts(X, 5)
#' chisq.test(tab)
split_counts <- function (samples, k = 2) 
{
  if (k <= 0) 
    stop("k needs to be larger than 0")
  if (k >= length(samples)) 
    stop("k needs to be smaller than the number of samples")
  broken.up <- split(samples, cut(seq_along(samples), k, labels = FALSE))
  get_freq <- function(x){
    out <- data.frame(
      "0" = sum(x=="0"),
      "1" = sum(x=="1")
    )
    return(out)
  }
  raw.ns <- lapply(broken.up, get_freq)
  ns <- matrix(t(do.call(rbind, raw.ns)), ncol = k)
  colnames(ns) <- paste0("w=", 1:k)
  rownames(ns) <- c("0", "1")
  return(ns)
}
#' Row-normalise a matrix with non-negative entries
#'
#' @param mat a square matrix with non-negative entries.
#'
#' @return a row-stochastic matrix
#' @export norm_matrix
#'
#' @examples
#' M <- matrix(c(124, 3, 16,
#'         6, 109, 14,
#'         22, 9, 142),
#'       ncol = 3, nrow = 3, byrow = TRUE)
#'norm_matrix(M)
norm_matrix <- function(mat){
  t(apply(mat, 1, function(x) x/sum(x)))
}
#'  Compute the likelihood ratio statistic comparing a stationary Markov chain model for \code{seq} to a nonstationary model
#'
#' @param seq a vector containing observations from an stochastic process.
#' @param L an integer specifying how many windows will be compared
#' @param states the possible states of the Markov chain. If \code{NULL},
#' will be deduced from 'seq'.
#' @param eps a small number to avoid zero probabilities
#'
#' @return a likelihood ratio statistic comparing H1 = not stationary vs
#'  H0 = stationary.
#' @export compute_Hsquared
#' @examples
#' X <- rbinom(n = 1E5, size = 1, prob = 0.8)
#' compute_Hsquared(X)
compute_Hsquared <- function (seq, L = 10, states = NULL, eps = 1E-17) 
{
  ## page 190 of Thomas (2013)
  ##  "Some Simplifying Conditions for Markov Chain Modeling"
  if(is.null(states)) states <- sort(unique(seq))
  nstates <- length(states)
  csize <- length(seq)/L
  windows <- split(seq, ceiling(seq_along(seq)/csize))
  dt <- ps <- array(NA, dim = c(nstates, nstates, L))
  for(t in 1:L){
    dt[,, t] <- get_contingency_array(windows[[t]], 1, states = states)
    ps[,, t] <- norm_matrix(dt[,, t])
  }
  pooledCounts <- apply(dt, c(1, 2), sum)
  pooledP <- norm_matrix(pooledCounts)
  Lsep <- 0
  Lpool <- 0
  for(t in 1:L){
    for(i in 1:nstates){
      for(j in 1:nstates){
        Lsep <- Lsep + dt[i, j, t] * log(ps[i, j, t] + eps)
        Lpool <- Lpool +  dt[i, j, t] * log(pooledP[i, j] + eps)
      }
    }
  }
  out <- 2*(Lsep - Lpool)
  return(out)
}
#' Compute the Bayes factor testing stationarity of a Markov chain.
#'
#' @param seq a vector containing observations from an stochastic process.
#' @param L an integer specifying how many windows will be compared
#' @param states the possible states of the Markov chain. If \code{NULL},
#' will be deduced from 'seq'.
#'
#' @return the approximate Bayes factor testing H1 = not stationary vs
#'  H0 = stationary.
#' @export compute_stationarity_BF
#'
#' @examples
#' X <- rbinom(n = 1E5, size = 1, prob = 0.8)
#' compute_stationarity_BF(X)
compute_stationarity_BF <- function(seq, L = 10, states = NULL){
  n <- length(seq)
  if(is.null(states)) states <- sort(unique(seq))
  Hsq <- compute_Hsquared(seq = seq, L = L, states = states)
  BF <- (Hsq - 2*log(n))/2
  return(BF)
}