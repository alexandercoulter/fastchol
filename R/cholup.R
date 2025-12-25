#' Cholesky Rank-1 Update
#'
#' @description
#' This function calculates the rank-1 update of Cholesky factor `CF`, e.g.
#' the Cholesky factor of matrix `(CF)(CF') + xx'`. Here, `CF` is a square `p`
#' \eqn{\times} `p` Cholesky factor, and `x` is any arbitrary vector of length
#' `p`. The user may specify whether `CF` is lower-triangular (in which case
#' `lower = TRUE`) or upper-triangular (in which case `lower = FALSE`). This
#' function applies Givens rotations to triangularize the joint matrix `(CF x)`;
#' due to cache latency, it will generally be faster when `CF` is the lower
#' Cholesky factor.
#' 
#' @param CF Cholesky factor of dimension `p` \eqn{\times} `p`
#' @param x Any `p`-long vector
#' @param lower Boolean, dictating whether `CF` is the lower Cholesky factor
#' (`TRUE`, default) or upper Cholesky factor (`FASLE`). Also dictates whether
#' the returned Choleskyfactor is lower- or upper-triangular.
#'
#' @returns Cholesky factor of (e.g.) `(CF)(CF') + xx'`, where `(CF)(CF')` is
#' the positive definite matrix associated with `CF`.
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
#' # Generate positive definite matrix
#' p = 20
#' M = diag(p) + tcrossprod(matrix(rnorm(p * 2), p, 2))
#' 
#' # Generate vector for rank-1 update
#' x = rnorm(p)
#' 
#' # Calculate Cholesky factor of M
#' U = chol(M)
#' 
#' # Calculate rank-1 update
#' Ux = cholup(U, x, lower = FALSE)
#' 
#' # Check against Cholesky factor of U'U + xx'
#' Mx = M + tcrossprod(x)
#' max(abs(Ux - chol(Mx)))
cholup = function(CF, x, lower = TRUE){
  
  p = length(x)
  R = CF
  R[1] = R[1] + 0
  x0 = x
  x0[1] = x0[1] + 0
  if(lower) cholupL_Rcpp(R, x0) else cholupU_Rcpp(R, x0)
  R
  
}