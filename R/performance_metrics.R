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
#' @seealso  \code{\link[BinaryMarkovChains]{split_switching_ratio}} 
#' \code{\link[BinaryMarkovChains]{MC_neff}} \code{\link[BinaryMarkovChains]{scaled_alpha}}
#'  \code{\link[BinaryMarkovChains]{scaled_average_transitions}}
#'  @examples 
#' X1 <- c(0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
#' X2 <- c(0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
#' switching_ratio(X1)
#' switching_ratio(X2)
#' switching_ratio(X1[1:8])
#' switching_ratio(X1[9:16])
#' switching_ratio(X2[1:8])
#' switching_ratio(X2[9:16])
switching_ratio <- function(samples){
  N <- length(samples)
  X.0 <- as.character(samples[1])
  occ.time <- sum(samples)
  max.transitions <- get_maxK(S = occ.time, M = N, X0 = X.0)
  Sy <- sum(form_Y(samples)) 
  delta <- Sy/max.transitions
  if(max.transitions <= 0) delta <- 0
  return(delta)
}

#' Split-Scaled number of state changes
#'
#' Computes the ratio between the number of observed state changes and
#'  the maximum possible number of state changes given the observed occupation time in 
#'  a number \code{k} of windows.
#' @param samples a collection of N samples from a two-state Markov chain.
#' @param  k a number of windows into which to break up the samples.
#' Default is \code{k=2}.
#'
#' @return a list of scaled number of state transitions (between 0 and 1) in each of the \code{k} windows.
#' @export split_switching_ratio
#' @seealso  \code{\link[BinaryMarkovChains]{switching_ratio}} \code{\link[BinaryMarkovChains]{MC_neff}}
#'  \code{\link[BinaryMarkovChains]{scaled_alpha}} \code{\link[BinaryMarkovChains]{scaled_average_transitions}}
#' @details This function breaks up the samples into \code{k} windows computes \code{switching_ratio} for
#' each one. Works in the same spirit as split-Rhat and other
#'  metrics which compare one part of the chain to the others.
#' @examples 
#' X <- rbinom(1E6, size = 1, prob = .2)
#' split_switching_ratio(X, k = 2)
#' split_switching_ratio(X, k = 11)
#' plot(split_switching_ratio(X, k = 200), type = "l")
split_switching_ratio <- function(samples, k = 2){
  if(k <= 0) stop("k needs to be larger than 0")
  if(k >= length(samples)) stop("k needs to be smaller than the number of samples")
  broken.up <- split(samples, cut(seq_along(samples), k, labels = FALSE))
  raw.deltas <- lapply(broken.up, switching_ratio)
  deltas <- unlist(raw.deltas)
  ns <- unlist(lapply(broken.up, length))
  names(deltas) <- paste0("n=", ns)
  return(deltas)
}