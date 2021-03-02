#' Autocorrelation function for a two-state Markov chain.
#'
#' @param lag the lag at which to compute autocorrelation. Should normally be non-negative integer.
#' @param alpha a transition probability (between 0 and 1).
#' @param beta a transition probability (between 0 and 1).
#'
#' @return theoretical autocorrelation at \code{lag}
#' @export
#'
#' @examples
#' library(markovchain)
#' true.alpha <- .02
#' true.beta <- .2
#' mat <- abs(matrix(c(1-true.alpha, true.alpha, true.beta, 1-true.beta),
#'  ncol = 2, nrow = 2, byrow = TRUE))
#' MC.binary <- new("markovchain", states = c("0", "1"),
#'  transitionMatrix = mat, name = "Binary")
#' M <- 1E4
#' outs <- markovchainSequence(n = M, markovchain = MC.binary, t0 = "0")
#' X <- as.numeric(outs)
#' ## Comparing estimated and theoretical autocorrelations
#' spec <- acf(X, plot = FALSE)
#' true.rhos <- MC_autocorr(lag = 0:(length(spec$acf)-1), alpha = true.alpha, beta = true.beta)
#' plot(true.rhos, as.numeric(spec$acf))
#' abline(a = 0, b = 1, lwd = 2, lty = 2)
MC_autocorr <- function(lag, alpha, beta){
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1")
  if(beta < 0 || beta > 1) stop("beta must be between 0 and 1")
  ab <- alpha + beta
  if(ab < 1){
    ans <- exp(lag * log1p(-ab))
  }else{
    ans <- (1-ab)^lag
  }

 return(ans)
}
MC_autocorr <- Vectorize(MC_autocorr)
