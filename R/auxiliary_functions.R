#' Vector of state-transitions.
#'
#' Constructs the binary vector of state transitions given a vector \code{x}.
#' @param x vector of length M containing realisations of a two-state (binary) discrete-time Markov chain.
#'
#' @return vector of state-transitions (changes).
#' @export form_Y
#'
#' @examples
#' X <- c(0, 1, 0, 1, 1, 0, 1)
#' Y <- form_Y(X)
#' sum(Y) # number of state-transitions (changes)
form_Y <- function(x){
  return(
    c(0, abs(diff(x)))
  )
}
#
#' Maximum alpha for a given marginal equilibrium probability p.
#'
#' @param p a probability of success in a two-state Markov chain.
#'
#' @return maximum alpha given p, alpha_m.
#' @export get_max_alpha
#'
#' @examples
#' get_max_alpha(.1)
#' get_max_alpha(.5)
#' get_max_alpha(.9)
get_max_alpha <- function(p){
  if(p < 0 || p > 1) stop("p must be between 0 and 1")
  ans <- min(p/(1-p), 1)
  return(ans)
}
#' Compute transition probability beta from alpha and p
#'
#' @param alpha 
#' @param p 
#'
#' @return beta
#' @export get_beta
#'
#' @examples
#' p <- 0.2
#' alpha <- get_max_alpha(p)/2
#' beta <- get_beta(alpha = alpha, p = p)
get_beta <- function(alpha, p){
  beta <- exp(log(alpha) + log1p(-p) - log(p))
  return(beta)
}
#' The a parameter in the state-change Markov chain.
#'
#' @param alpha a transition probability (between 0 and 1).
#' @param beta a transition probability (between 0 and 1).
#'
#' @return a transition probability.
#' @importFrom matrixStats logSumExp
#' @export get_a
#' @details  If X_k is a two-state Markov chain with transition probabilities \code{alpha} and \code{beta}, then
#'  the state-changes Y_k = |X_k - X_{k-1}| also form a two-state Markov chain with transition probabiliti`es
#'  \code{a} and \code{b}.
#'
get_a <- function(alpha, beta){
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1")
  if(beta < 0 || beta > 1) stop("beta must be between 0 and 1")
  p <- alpha/(alpha + beta)
  lnum <- matrixStats::logSumExp(c(
    log(alpha) + log1p(-alpha) + log1p(-p),
    log(beta) + log1p(-beta) + log(p)
  ))
  ldenom <- matrixStats::logSumExp(c(
    log1p(-beta) + log(p),
    log1p(-p) + log1p(-alpha)
  ))
  ans <- exp(lnum-ldenom)
  return(ans)
}
#' The b parameter in the state-change Markov chain.
#'
#' @param alpha a transition probability (between 0 and 1).
#' @param beta a transition probability (between 0 and 1).
#'
#' @return b transition probability.
#' @importFrom matrixStats logSumExp
#' @export get_b
#' @details  If X_k is a two-state Markov chain with transition probabilities \code{alpha} and \code{beta}, then
#'  the state-changes Y_k = |X_k - X_{k-1}| also form a two-state Markov chain with transition probabiliti`es
#'  \code{a} and \code{b}.
#'
get_b <- function(alpha, beta){
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1")
  if(beta < 0 || beta > 1) stop("beta must be between 0 and 1")
  p <- alpha/(alpha + beta)
  lnum <- matrixStats::logSumExp(c(
    log1p(-beta) + log(alpha) + log1p(-p),
    log1p(-alpha) + log(beta) + log(p)
  ))
  ldenom <- matrixStats::logSumExp(c(
    log(beta) + log(p),
    log(alpha) + log1p(-p)
    ))
  ans <- exp(lnum-ldenom)
  return(ans)
}
#'  Maximum number of transitions given occupation time, total iterations and initial state.
#'
#' @param Sx the sum of 1's in M samples.
#' @param M number of samples.
#' @param X0 initial state. Either \code{X0 == "0"} or \code{X0 == "1"}.
#'
#' @return maximum possible number of state-transitions (changes).
#' @export get_maxK
#'
#' @examples
#' get_maxK(Sx = 3, M = 10) # should be 6
get_maxK <- function(Sx, M, X0 = "0"){
  m <- min(Sx, M-Sx)
  if(X0 == "0"){
    if(Sx == M){
      warning("Sum is M but X0 = 0. Something is not right.")
      maxT <- 0
    }else{
      maxT <- ifelse(m == Sx, 2*m, 2*m -1)
    }
  }else{
    if(Sx == 0){
      warning("Sum is 0 but X0 = 1. Something is not right.")
      maxT <- 0
    }else{
      maxT <- ifelse(m == Sx, 2*m-1, 2*m)
    }
  }
  return(maxT)
}
get_maxK <- Vectorize(get_maxK)
