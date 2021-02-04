#' Maximum _a posteriori_ estimate of the alpha transition probability
#'
#' @param dmat a 2 x 2 matrix with observed transitions.
#' @param k shape1 parameter of Beta prior on alpha (k > 0). Default is \code{k = 1}.
#' @param v shape2 parameter of Beta prior on alpha (v > 0). Default is \code{v = 1}.
#' @param p marginal success probability.
#'
#' @return MAP estimate of alpha
#' @export get_alpha_map
#'
#' @examples
#' library(markovchain)
#' X <- sample(c(0, 1), 10000, replace = TRUE)
#' obs.trans <- createSequenceMatrix(X, sanitize = FALSE, possibleStates = c("0", "1"))
#' get_alpha_map(dmat = obs.trans, p = 1/2) ## should be close to 1/2 because of independent sampling
get_alpha_map <- function(dmat, k = 1, v = 1, p){
  z <- (1-p)/p
  a <- dmat[1, 1] + v - 1
  b <- dmat[1, 2] + dmat[2, 1]
  c <- dmat[2, 2] + k - 1
  raw <- -(sqrt((c^2+2*b*c+b^2)*z^2+((2*a-2*b)*c-2*b^2-2*a*b)*z+b^2+2*a*b+a^2)+(-c-b)*z-b-a)/((2*c+2*b+2*a)*z)
  if(is.nan(raw)) raw <- 0
  ans <- min(max(0, raw), 1)
  return(ans)
}
