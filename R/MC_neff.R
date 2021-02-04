#' Compute the theoretical effective sample size in two-state Markov chain.
#'
#'This is the theoretical effective sample size for a sample of size \code{N} given parameters \code{alpha} and \code{p}.
#' @param N number of iterations
#' @param alpha a transition probability (between 0 and 1)
#' @param p (marginal) success probability
#'
#' @return effective sample size
#' @export MC_neff_theoretical
#' @seealso  \code{\link[BinaryMarkovChains]{MC_neff}}
#' @examples
#' p <- .25
#' MC_neff_theoretical(N = 100, alpha = p, p = p) ## independent sampling gives N_eff = N
MC_neff_theoretical <- function(N, alpha, p){
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1")
  if(p < 0 || p > 1) stop("p must be between 0 and 1")
  sr <- 2*(p-alpha)/alpha
  if(sr < -1){
    ans <- NaN
  }else{
    ans <- exp(log(N) - log1p(sr))
  }
  return(ans)
}
#' Estimate of the effective sample for two-state Markov chains.
#'
#' @param samples a vector of size N containing realisations of a two-state Markov chain.
#' @param p (optional) marginal success probability. If \code{p = NULL}, the sample mean will be used.
#'
#' @return estimated effective sample size
#' @importFrom markovchain createSequenceMatrix
#' @export MC_neff
#' @seealso  \code{\link[BinaryMarkovChains]{MC_neff_theoretical}}
#' \code{\link[BinaryMarkovChains]{scaled_alpha}}
#' \code{\link[BinaryMarkovChains]{scaled_average_transitions}}
#' \code{\link[BinaryMarkovChains]{switching_ratio}}
#' @details  This functions provides an estimate of the effective sample size by plugging in an estimate of alpha in the form of its
#' maximum _a posteriori_  estimate (see \code{\link[BinaryMarkovChains]{get_alpha_map}}).
#' If p is not known (\code{p = NULL}) then the sample mean, \code{mean(samples)} is used as an estimate.
#' @examples
#' X <- sample(c(0, 1), 1000, replace = TRUE)
#' MC_neff(samples = X, p = 1/2)
#' MC_neff(samples = X)
MC_neff <- function(samples, p = NULL){
  if(is.null(p)) p <- mean(samples, na.rm = TRUE)
  if(p < 0 || p > 1) stop("p must be between 0 and 1")
  obs.transitions <- markovchain::createSequenceMatrix(samples, sanitize = FALSE,
                                                       possibleStates = c("0", "1"))
  alpha.hat <- get_alpha_map(dmat = obs.transitions, k = 1, v = 1, p = p)
  ans <- MC_neff_theoretical(N = length(samples), alpha = alpha.hat, p = p)
  return(ans)
}
