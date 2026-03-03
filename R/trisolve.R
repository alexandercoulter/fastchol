#' Triangular Solver
#'
#' @description
#' This function is a custom implementation of the base `R` functions 
#' `forwardsolve` (when `lower = TRUE`) and `backsolve` (when `lower = FALSE`),
#' evaluating \eqn{(CF)^{-1}X}.
#'
#' @param CF Triangular matrix of dimension `n` \eqn{\times} `n`
#' @param X Any matrix of dimension `n` \eqn{\times} `p`, or vector of length `n`
#' @param lower Boolean, dictating whether `CF` is lower-triangular
#' (`TRUE`, default) or upper-triangular (`FASLE`).
#'
#' @returns Solution for `Y` in linear system `(CF)Y = X`, i.e. \eqn{(CF)^{-1}X}
#' @export
#'
#' @examples
#' # Generate positive definite matrix
#' n = 20
#' M = diag(n) + tcrossprod(matrix(rnorm(n * 2), n, 2))
#' 
#' # Calculate Cholesky factors
#' U = chol(M)
#' L = t(U)
#' 
#' # Generate RHS matrix for linear system
#' p = 5
#' X = matrix(rnorm(n * p), n, p)
#' 
#' # Calculate triangular solves
#' Lx = trisolve(L, X, lower = TRUE)
#' Ux = trisolve(U, X, lower = FALSE)
#' 
#' # Compare to forwardsolve and backsolve
#' max(abs(Lx - forwardsolve(L, X)))
#' max(abs(Ux - backsolve(U, X)))
trisolve = function(CF, X, lower = TRUE){
  
  X0 = X
  X0[1] = X0[1] + 0
  
  if(is.matrix(X0)){
    
    if(lower) LsolveX_Rcpp(CF, X0) else UsolveX_Rcpp(CF, X0)
    
  } else {
    
    if(lower) Lsolvex_Rcpp(CF, X0) else Usolvex_Rcpp(CF, X0)
    
  }
  
  X0
  
}