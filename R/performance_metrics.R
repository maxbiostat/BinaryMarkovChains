#' Scaled average state changes.
#'
#' @param samples a collection of N samples from a two-state Markov chain.
#' @param p (optional) marginal success probability.
#'
#' @return normalised number of state transitions (changes).
#' @export scaled_average_transitions
#' @details The function computes the mean number of state changes, \code{y.hat = Sy/N}, where \code{Sy} is
#' the sum of state changes and \code{N} is the total number of iterations.
#' The function then returns  \code{y.hat} scaled by its maximum theoretical value given \code{p}.
#' If p is not known (\code{p = NULL}) then the sample mean, \code{mean(samples)} is used as an estimate.
#' @seealso  \code{\link[BinaryMarkovChains]{MC_neff}} \code{\link[BinaryMarkovChains]{scaled_alpha}} \code{\link[BinaryMarkovChains]{switching_ratio}}
scaled_average_transitions <- function(samples, p = NULL){
   if(is.null(p)) p <- mean(samples, na.rm = TRUE)
   y.hat <- mean(form_Y(samples), na.rm = TRUE)
   alpha.max <- get_max_alpha(p)
   ld <- log(2) + log(alpha.max) + log1p(-p)
   phi <- exp(log(y.hat)-ld)
   return(phi)
 }

#' Scaled alpha estimate.
#'
#' @param samples a collection of N samples from a two-state Markov chain.
#' @param k shape1 parameter of Beta prior on alpha (k > 0). Default is \code{k = 1}.
#' @param v shape2 parameter of Beta prior on alpha (v > 0). Default is \code{v = 1}.
#' @param p (optional) marginal success probability.
#'
#' @return estimate of alpha scaled by its maximum value given \code{p}.
#' @importFrom  markovchain createSequenceMatrix
#' @export scaled_alpha
#' @details Function computes the maximum _a posteriori_ estimate of the transition probability alpha and then
#' scales that estimate by its maximum theoretical value given \code{p}.
#' If p is not known (\code{p = NULL}) then the sample mean, \code{mean(samples)} is used as an estimate.
#' @seealso  \code{\link[BinaryMarkovChains]{MC_neff}}
#'  \code{\link[BinaryMarkovChains]{scaled_average_transitions}}
#'  \code{\link[BinaryMarkovChains]{switching_ratio}}
scaled_alpha <- function(samples, k = 1, v = 1, p = NULL){
  if(is.null(p)) p <- mean(samples, na.rm = TRUE)
  alpha.max <- get_max_alpha(p)
  data.mat <- markovchain::createSequenceMatrix(samples, sanitize = FALSE,
                                                possibleStates = c("0", "1"))
  alpha.hat <- get_alpha_map(dmat = data.mat, k = k, v = v, p = p)
  psi <- exp(log(alpha.hat)-log(alpha.max))
 return(psi)
}

#' Scaled number of state changes
#'
#' Computes the ratio between the number of observed state changes and
#'  the maximum possible number of state changes given the observed occupation time.
#' @param samples a collection of N samples from a two-state Markov chain.
#'
#' @return scaled number of state transitions (between 0 and 1).
#' @export switching_ratio
#' @details This index is similar to \code{\link[BinaryMarkovChains]{scaled_average_transitions}} but
#'  does not require any knowledge about any of the Markov chain parameters.
#'  On the other hand, it is computed **conditional** on the observed occupation time and thus inherits its uncertainty.
#' @seealso  \code{\link[BinaryMarkovChains]{MC_neff}} \code{\link[BinaryMarkovChains]{scaled_alpha}} \code{\link[BinaryMarkovChains]{scaled_average_transitions}}
switching_ratio <- function(samples){
  N <- length(samples)
  occ.time <- sum(samples)
  max.transitions <- get_maxK(S = occ.time, M = N)
  Sy <- sum(form_Y(samples))
  delta <- Sy/max.transitions
  return(delta)
}
