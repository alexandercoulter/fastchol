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
#'
#' @returns Cholesky factor of (e.g.) `(CF)(CF') + xx'`, where `(CF)(CF')` is
#' the positive definite matrix associated with `CF`.
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
cholupL = function(CF, x, lower = TRUE){
  
  p = length(x)
  R = CF
  R[1] = R[1] + 0
  x0 = x
  x0[1] = x0[1] + 0
  if(lower) cholupL_Rcpp(R, x0) else cholupU_Rcpp(R, x0)
  R
  
}