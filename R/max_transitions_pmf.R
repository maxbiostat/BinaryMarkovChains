#' Probability mass function for the maximum number of transitions.
#'
#' Computes the p.m.f. of the maximum possible number of state transitions (Delta_m) in a
#'  two-state Markov chain **conditioning** on X_0 = 0 or X_0 = 1.
#' @param x maximum number of state transitions in \code{n} iterations.
#' @param n number of Markov chain iterations.
#' @param alpha a transition probability (between 0 and 1).
#' @param beta a transition probability (between 0 and 1).
#' @param X0 initial state. Either \code{X0 == "0"} or \code{X0 == "1"}.
#' @param log logical if TRUE, probabilities p are given as log(p).
#'
#' @return  (log) p, where p = Pr(Delta_m = x)
#' @export max_transitions_pmf
#' @seealso  \code{\link[BinaryMarkovChains]{occupation_time_pmf}}, \code{\link[BinaryMarkovChains]{recurrence_time_pmf}}
#' @examples
#' M <- 100
#' deltas.0 <- get_maxK(Sx = 0:(M-1), M = M, X0 = "0")[1:M]
#' deltas.1 <- get_maxK(Sx = 1:M, M = M, X0 = "1")[1:M]
#' probs.X0 <- max_transitions_pmf(x = deltas.0, n = M, alpha = .02, beta = .2, X0 = "0")
#' probs.X1 <- max_transitions_pmf(x = deltas.1, n = M, alpha = .02, beta = .2, X0 = "1")
#' sum(probs.X0); sum(probs.X1) # should both be 1
#' plot(deltas.0, probs.X0, type = "h", lwd = 3, ylab = "Probability", xlab = expression(Delta[m]), main = "Maximum  transitions")
#' points(deltas.1, probs.X1, type = "h", lwd = 3, col = 2)
#' legend(x = "topright", col = 1:2, legend = c(expression(X[0]==0),
#'  expression(X[0]==1)), bty = 'n', pch = 15)
max_transitions_pmf <- function(x, n , alpha, beta, X0 = c("0", "1"), log = FALSE){
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1")
  if(beta < 0 || beta > 1) stop("beta must be between 0 and 1")
  if(x > n) stop("x needs to be less than or equal to n")
  l1 <- occupation_time_pmf(x = ceiling(x/2), n = n, alpha = alpha, beta = beta, X0 = X0, log = TRUE)
  l2 <- occupation_time_pmf(x = ceiling(((2*n -1)-x)/2), n = n, alpha = alpha, beta = beta, X0 = X0, log = TRUE)
  if(x == 0){
    ans <- matrixStats::logSumExp(c(l1, l2))
  }else{
    ans <- matrixStats::logSumExp(c(l1, l2)) -log(2)
  }
  if(!log) ans <- exp(ans)
  return(ans)
}
max_transitions_pmf <- Vectorize(max_transitions_pmf)
