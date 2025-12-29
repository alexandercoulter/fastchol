#' Cholesky Rank-K Update
#'
#' @description
#' This function calculates the rank-K update of Cholesky factor `CF`, e.g.
#' the Cholesky factor of matrix `(CF)(CF') + XX'`. Here, `CF` is a square `p`
#' \eqn{\times} `p` Cholesky factor, and `X` is any arbitrary matrix with `p`
#' rows and `K` columns The user may specify whether `CF` is lower-triangular
#' (in which case `lower = TRUE`) or upper-triangular (in which case
#' `lower = FALSE`). This function applies Givens rotations to triangularize the
#' joint matrix `(CF x)`; due to cache latency, it will generally be faster when
#' `CF` is the lower Cholesky factor.
#' 
#' @param CF Cholesky factor of dimension `p` \eqn{\times} `p`
#' @param X Any `p` \eqn{\times} `K` matrix
#' @param lower Boolean, dictating whether `CF` is the lower Cholesky factor
#' (`TRUE`, default) or upper Cholesky factor (`FASLE`). Also dictates whether
#' the returned Choleskyfactor is lower- or upper-triangular.
#'
#' @returns Cholesky factor of (e.g.) `(CF)(CF') + xx'`, where `(CF)(CF')` is
#' the positive definite matrix associated with `CF`.
#' @export
#'
#' @examples
#' # Generate positive definite matrix
#' p = 20
#' M = diag(p) + tcrossprod(matrix(rnorm(p * 2), p, 2))
#' 
#' # Generate matrix for rank-K update, i.e. M + XX' = M + sum_k {X[ , k] * X[ , k]'}
#' K = 3
#' X = matrix(rnorm(p * K), p, K)
#' 
#' # Calculate Cholesky factor of M
#' U = chol(M)
#' 
#' # Calculate rank-1 update
#' UX = cholupK(U, X, lower = FALSE)
#' 
#' # Check against Cholesky factor of U'U + XX'
#' MX = M + tcrossprod(X)
#' max(abs(UX - chol(MX)))
cholupK = function(CF, X, lower = TRUE){
  
  R = CF
  R[1] = R[1] + 0
  X0 = X
  X0[1] = X0[1] + 0
  if(lower) cholupKL_Rcpp(R, X0) else cholupKU_Rcpp(R, X0)
  R
  
}