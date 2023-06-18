#' Joint probability between two indicators
#'
#' @param a a vector with an indicator variable.
#' @param b a vector with an indicator variable.
#' @param log logical. If \code{TRUE} returns the log probability. 
#'
#' @return joint probability
#' @export joint_p
#'
#' @examples
#' X <- rbinom(1000, size = 1, p = .25)
#' Y <- rbinom(1000, size = 1, p = .75)
#' joint_p(X, Y)
joint_p <- function(a, b, log = FALSE){
  if(length(a) != length(b)) stop("Vectors not of the same size")
  ans <- log(mean(a*b))
  if(!log) ans <- exp(ans)
  return(ans)
}
#' Covariance between binary variables.
#'
#' @param pij joint probability.
#' @param pi marginal probability of the first variable.
#' @param pj marginal probability of the second variable.
#'
#' @return covariance.
#' @export bin_cov
#'
#' @examples
#' bin_cov(pij = 0.5, pi = .2, pj = .2)
#' bin_cov(pij = 0.0, pi = .2, pj = .2)
bin_cov <- function(pij, pi, pj){
  return( pij - exp(log(pi) + log(pj)) )
}
#' Correlation between binary variables#'
#' @param pij joint probability.
#' @param pi marginal probability of the first variable.
#' @param pj marginal probability of the second variable.
#'
#' @return correlation coefficient.
#' @export bin_corr
#'
#' @examples
#' bin_corr(pij = 0.1, pi = 0.1, pj = .2) 
#' bin_corr(pij = 0.0, pi = .2, pj = .2)
#' # bin_corr(pij = 0.5, pi = .2, pj = .2) ## Error
bin_corr <- function(pij, pi, pj){
  num <- pij - exp(log(pi) + log(pj))
  sgn <- sign(num)
  lnum <- log(abs(num))
  ldenom <- 0.5 * ( log(pi) + log1p(-pi) + log(pj) + log1p(-pj) )
  ans <- sgn * exp(lnum - ldenom)
  if(ans > 1)  stop("Computed correlation is larger than 1, check your probabilities")
  return(ans)
}
#' Covariance matrix for binary variables
#'
#' @param x a matrix or data.frame containing only binary variables.
#' @param ncores (optional) number of cores to be used.
#'
#' @return covariance matrix.
#' @importFrom stats var
#' @export binary_cov_matrix
#'
#' @examples
#' N <- 100
#' M <- 1000
#' X <- matrix(rbinom(n = N*M, size = 1, p = .234), ncol = N, nrow = M)
#' covmat.binary <- binary_cov_matrix(X)
#' covmat.default <- var(X)
#' plot(covmat.binary ~ covmat.default)
#' abline(a = 0, b = 1, lty = 2)
binary_cov_matrix <- function(x, ncores = 2){
  Var1 <- Var2 <- NULL
  Vs <- apply(x, 2, var)
  if(any(Vs ==0)) stop("There are zero-variance variables!")
  N <- ncol(x)
  grid <- subset(expand.grid(1:N, 1:N), Var1 < Var2)
  ps <- as.numeric(colMeans(x))
  covs <- unlist(
    parallel::mclapply(1:nrow(grid), function(row) {
      i <- grid[row, 1]
      j <- grid[row, 2]
      ans <- bin_cov(pij = joint_p(x[, i], x[, j]), pi = ps[i], pj = ps[j])
      return(ans)
    }, mc.cores = ncores)
  ) 
  all.comps <- cbind(grid, covs)
  cov.mat <- matrix(NA, ncol = N, nrow = N)
  diag(cov.mat) <- Vs
  for(k in 1:nrow(all.comps)){
    i <- all.comps[k, 1]
    j <- all.comps[k, 2]
    cov.mat[i, j] <- cov.mat[j, i] <- all.comps[k, 3]
  }
  return(cov.mat)
}
#' Correlation matrix for binary variables.
#'
#' @param x a matrix or data.frame containing only binary variables.
#' @param ncores (optional) number of cores to be used.
#'
#' @return a correlation matrix.
#' 
#' @importFrom stats var
#' @export binary_corr_matrix
#'
#' @examples
#' N <- 100
#' M <- 1000
#' X <- matrix(rbinom(n = N*M, size = 1, p = .234), ncol = N, nrow = M)
#' cormat.binary <- binary_corr_matrix(X)
#' cormat.default <- cor(X)
#' plot(cormat.binary ~ cormat.default)
#' abline(a = 0, b = 1, lty = 2)
binary_corr_matrix <- function(x, ncores = 2){
  Var1 <- Var2 <- NULL
  Vs <- apply(x, 2, var)
  if(any(Vs ==0)) stop("There are zero-variance variables!")
  N <- ncol(x)
  grid <- subset(expand.grid(1:N, 1:N), Var1 < Var2)
  ps <- as.numeric(colMeans(x))
  corrs <- unlist(
    parallel::mclapply(1:nrow(grid), function(row) {
      i <- grid[row, 1]
      j <- grid[row, 2]
      ans <- bin_corr(pij = joint_p(x[, i], x[, j]), pi = ps[i], pj = ps[j])
      return(ans)
    }, mc.cores = ncores)
  ) 
  all.comps <- cbind(grid, corrs)
  corr.mat <- matrix(NA, ncol = N, nrow = N)
  diag(corr.mat) <- 1
  for(k in 1:nrow(all.comps)){
    i <- all.comps[k, 1]
    j <- all.comps[k, 2]
    corr.mat[i, j] <- corr.mat[j, i] <- all.comps[k, 3]
  }
  return(corr.mat)
}

#' Multivariate sampling efficiency using the log-determinant
#'
#' @param varMat covariance matrix
#' @param covMat long term covariance matrix
#'
#' @return multivariate efficiency
#' @export m_eff_det
#' @details
#' The function uses \code{base::determinant()} to get the required determinants
#' 
#'
m_eff_det <- function(varMat, covMat) {
  p <- ncol(covMat)
  det.L <- determinant(varMat)
  det.S <- determinant(covMat)
  eff <- 
    exp((as.numeric(det.L$modulus) - as.numeric(det.S$modulus)) / p)
  return(eff)
}

#' Multivariate sampling efficiency using (log) eigen values
#'
#' @param varMat covariance matrix
#' @param covMat long term covariance matrix
#'
#' @return multivariate efficiency
#' @export m_eff_eigen
#' @details
#' The function uses \code{base::eigen()} to get the required determinants
#' we take the log of the eigenvalues, and sum them to get the product in log-space
#'
m_eff_eigen <- function(varMat, covMat) {
  p <- ncol(varMat)
  log.det.var.p <- sum(log(eigen(
    varMat, symmetric = TRUE,
    only.values = TRUE
  )$values))
  log.det.covmat.p <-
    sum(log(eigen(covMat, only.values = TRUE)$values))
  ess <- exp((log.det.var.p - log.det.covmat.p) / p)
  return(ess)
}

#' Multivariate ess
#'
#' @param N number of draws
#' @param Lambda covariance matrix
#' @param Sigma long-term covariance matrix
#' @param eigen boolean. If \code{TRUE}, the eigenvalue method is used. Otherwise,
#' the default \code{FALSE} uses the determinant (LU decomposition) method.
#'
#' @return the multivariate effective sample size
#' @export theo_multivariate_ESS
#'
theo_multivariate_ESS <- function(N, Lambda, Sigma, eigen = FALSE) {
  p <- ncol(Sigma)
  if (nrow(Sigma) != p)
    stop("Sigma is not square")
  if (nrow(Lambda) != ncol(Lambda))
    stop("Lambda is not square")
  if (nrow(Lambda) != p)
    stop("Lambda is not the same dimension as Sigma")
  
  if(eigen){
    eff <- m_eff_eigen(varMat = Lambda, covMat  = Sigma)
  }else{
    eff <- m_eff_det(varMat = Lambda, covMat  = Sigma)
  }
  return(N * eff)
}

#' Multivariate ESS for binary variables
#'
#' @param x a matrix or data.frame containing only binary variables
#' @param covmat (optional) covariance matrix.
#' @param ncores (optional) number of cores to be used.
#' @param ... arguments for mcse.multi function.
#'  Don't use this if a suitable matrix estimate from mcse.multi or
#'  mcse.initseq is already obtained (copied from **mcmcse** documentation).
#'
#' @return the multivariate ESS of the variables.
#' 
#' @importFrom  mcmcse mcse.multi
#' @export binary_multiESS
#' @details This is strongly based on the function \code{multiESS} of the **mcmcse** package.
#'
#' @examples
#' N <- 100
#' M <- 1000
#' X <- matrix(rbinom(n = N*M, size = 1, p = .234), ncol = N, nrow = M)
#' binary_multiESS(X)
#' mcmcse::multiESS(X)
binary_multiESS <- function (x, covmat = NULL,
                             ncores = 2, eigen = FALSE, ...) 
{
  chain <- as.matrix(x)
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("'x' must be a matrix or data frame.")
  n = dim(chain)[1]
  p = dim(chain)[2]
  if (is.matrix(covmat)) {
    var_mat <- binary_cov_matrix(chain, ncores = ncores)
    ess <- theo_multivariate_ESS(N = n,
                            Lambda = var_mat, Sigma = covmat,
                            eigen = eigen)
  }
  else {
    covmat <- mcmcse::mcse.multi(chain, ...)$cov
    var_mat <- binary_cov_matrix(chain, ncores = ncores)
    ess <- theo_multivariate_ESS(N = n,
                            Lambda = var_mat, Sigma = covmat,
                            eigen = eigen)
  }
  return(ess)
}