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
get_transitions_windows <- function(seq, nw = 2){
  if (nw <= 0) 
    stop("nw needs to be larger than 0")
  if (nw >= length(seq)) 
    stop("nw needs to be smaller than the number of samples")
  states <- sort(unique(seq))
  nstates <- length(seq)
  broken.up <- split(seq, cut(seq_along(seq), nw, labels = FALSE))
  raw <- lapply(broken.up, get_contingency_array, 1)
  out <- array(NA, dim = c(nstates, nstates, nw))
  for(k in 1:nw) out[,, k] <- raw[[k]]
  return(out)
}
compute_LRS <- function(seq, nw = 2){
  states <- sort(unique(seq))
  nstates <- length(states)
  dt <-  get_transitions_windows(seq, nw)
  stats <- c()
  for(k in 1:nw){
    for( i in 1:nstates){
      for(j in 1:nstates){
        nijk <-  dt[i , j, k] + 1E-17
        nijp <- sum(dt[i, j, ]) + 1E-17
        npjk <- sum(dt[, j, k]) + 1E-17
        npjp <- sum(dt[, j, ]) + 1E-17
        s <- nijk * log((nijk/nijp)/(npjk/npjp))
        stats <- c(stats, s)
      }
    }
  }
  out <- 2*sum(stats)
  return(out)
}
