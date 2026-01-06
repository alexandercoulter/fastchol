#' Triangular Solver with Matrix Argument
#'
#' @description
#' This function is a custom implementation of the default execution of 
#' `forwardsolve` (when `lower = TRUE`) or `backsolve` (when `lower = FALSE`),
#' evaluating \eqn{(CF)^{-1}x}.
#'
#' @param CF Triangular matrix of dimension `p` \eqn{\times} `p`
#' @param x Any vector of length `p`.
#' @param lower Boolean, dictating whether `CF` is lower-triangular
#' (`TRUE`, default) or upper-triangular (`FASLE`).
#'
#' @returns Solution for `y` in linear system `(CF)y = x`, i.e. \eqn{(CF)^{-1}x}
#' @export
#'
#' @examples
#' # Generate positive definite matrix
#' p = 20
#' M = diag(p) + tcrossprod(matrix(rnorm(p * 2), p, 2))
#' 
#' # Calculate Cholesky factors
#' U = chol(M)
#' L = t(U)
#' 
#' # Generate vector for linear system
#' x = rnorm(p)
#' 
#' # Calculate triangular solves
#' Lx = trisolve(L, x, lower = TRUE)
#' Ux = trisolve(U, x, lower = FALSE)
#' 
#' # Compare to forwardsolve and backsolve
#' max(abs(Lx - forwardsolve(L, x)))
#' max(abs(Ux - backsolve(U, x)))
trisolveX = function(CF, X, lower = TRUE){
  
  X0 = X
  X0[1] = X0[1] + 0
  if(lower) LsolveX_Rcpp(CF, X0) else UsolveX_Rcpp(CF, X0)
  X0
  
}