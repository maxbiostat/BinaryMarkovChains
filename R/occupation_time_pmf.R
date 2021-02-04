#' Probability mass function for the occupation time.
#'
#' Computes the pmf for the occupation time (sum of '1's), S, in a two-state Markov chain **conditioning** on X_0 = 0 or X_0 = 1.
#' This is the expression given Equation (1) of Bhattacharya & Gupta (1980), with a minor typo correction.
#' @param x the number of '1's in n iterations; non-negative integer smaller than or equal to n.
#' @param n the number of iterations in Markov chain.
#' @param alpha a transition probability (between 0 and 1).
#' @param beta a transition probability (between 0 and 1).
#' @param X0 initial state, either \code{X0 == "0"} or \code{X0 == "1"}.
#' @param log logical if TRUE, probabilities p are given as log(p).
#'
#' @return (log) p, where p = Pr(S = x).
#' @importFrom matrixStats logSumExp
#' @export occupation_time_pmf
#' @references Bhattacharya, S. K., & Gupta, A. K. (1980). Occupation times for two-state Markov chains. Discrete Applied Mathematics, 2(3), 249-250.
#' @seealso \code{\link[BinaryMarkovChains]{recurrence_time_pmf}}, \code{\link[BinaryMarkovChains]{max_transitions_pmf}}
#' @examples
#' M <- 100 # number of Markov chain iterations
#' probs.X0 <- occupation_time_pmf(x = 0:M, n = M, alpha = .02, beta = .2, X0 = "0")
#' probs.X1 <- occupation_time_pmf(x = 0:M, n = M, alpha = .02, beta = .2, X0 = "1")
#' sum(probs.X0) # should both be 1
#' sum(probs.X1)
#' plot(0:M, probs.X0, type = "h", lwd = 3, xlab = "Occupation time")
#' points(0:M, probs.X1, type = "h", lwd = 3, col = 2)
#' legend(x = "topright", col = 1:2,
#'    legend = c(expression(X[0]==0), expression(X[0]==1)), bty = 'n', pch = 15)
#'
occupation_time_pmf <- function(x, n , alpha, beta, X0 = c("0", "1"), log = FALSE){
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1")
  if(beta < 0 || beta > 1) stop("beta must be between 0 and 1")
  if(x > n) stop("x needs to be less than or equal to n")
  k <- n-x ## x is the number of successes, this pmf is for failures (conditional on X_0 =failure)
  lp11 <- log1p(-alpha)
  lp21 <- log(beta)
  lp12 <- log(alpha)
  lp22 <- log1p(-beta)
  if(X0 == "0"){
    if(k == 0){
      ans <- lp12 + (n-1)*lp22
    }else{
      if(k == n){
        ans <- n*lp11
      }else{
        lterm1 <- (n-2*k-1)*lp22 + k *(lp12 + lp21)
        s_term <- function(r){
          out <- lchoose(k, r) + r *(lp11 + lp22 - lp12 - lp21) +
            matrixStats::logSumExp(c(lp22 + lchoose(n-k-1, k-r-1), lp12 + lchoose(n-k-1, k-r)))
          return(out)
        }
        s_term <- Vectorize(s_term)
        lterm2 <- matrixStats::logSumExp(s_term(0:k))
        ans <- lterm1 + lterm2
      }
    }
  }else{
    if(k == 0){
      ans <- n* lp22
    }else{
      if(k == n){
        ans <- lp21 + (n-1)*lp11
      }else{
        lterm1 <- (n-2*k)*lp22 + (k-1)*lp12 + k *lp21
        s_term_2 <- function(r){
          out <- lchoose(k-1, r) + r*(lp11 + lp22 - lp12 - lp21) +
            matrixStats::logSumExp(c(lp22 + lchoose(n-k, k-r-1), lp12 + lchoose(n-k, k-r)))
          return(out)
        }
        s_term_2 <- Vectorize(s_term_2)
        lterm2 <- matrixStats::logSumExp(s_term_2(0:(k-1)))
        ans <- lterm1 + lterm2
      }
    }
  }
  if(!log) ans <- exp(ans)
  return(ans)
}
occupation_time_pmf <- Vectorize(occupation_time_pmf)
