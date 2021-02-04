#' Probability mass function of the recurrence time.
#'
#' Computes the p.m.f of the recurrence time T  i.e., number of iterations to come
#' back to state 0 in a two-state Markov chain if on X_0 = 0.
#' Result is given in page 83 of Cox & Miller (1965).
#'
#' @param x number greater than 0 of iterations until hitting 0 again.
#' @param alpha a transition probability (between 0 and 1).
#' @param beta a transition probability (between 0 and 1).
#' @param log logical if TRUE, probabilities p are given as log(p).
#'
#' @return (log) p, where p = Pr(T = x).
#' @export recurrence_time_pmf
#' @seealso \code{\link[BinaryMarkovChains]{occupation_time_pmf}}, \code{\link[BinaryMarkovChains]{max_transitions_pmf}}
#' @references  Cox, D. R., & Miller, H. D. (1977). The Theory of Stochastic Processes (Vol. 134). CRC press.
#' @examples
#' times <- 1:50
#' timeprobs <- recurrence_time_pmf(x = times, alpha = .02, beta = .2)
#' plot(times, timeprobs, type = "h", lwd = 3,
#' ylab = "Probability", xlab = expression(T), main = "Recurrence time")
recurrence_time_pmf <- function(x, alpha, beta, log = FALSE){
  if(x <= 0) stop("x needs to be a positive integer")
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1")
  if(beta < 0 || beta > 1) stop("beta must be between 0 and 1")
  ans <- log(alpha) + log(beta) + (x-2)*log1p(-beta)
  if(!log) ans <- exp(ans)
  return(ans)
}
recurrence_time_pmf <- Vectorize(recurrence_time_pmf)
