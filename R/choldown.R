#' Cholesky Rank-1 Downdate
#'
#' @description
#' This function calculates the rank-1 downdate of lower/upper Cholesky factor
#' `CF`, i.e. the lower Cholesky factor of matrix `(CF)(CF') - xx'` (if e.g.
#' `CF` is lower-triangular), provided the resulting matrix is positive
#' definite. This function applies hyperbolic rotations to triangularize the
#' joint matrix `(CF x)`. The user may specify whether `CF` is lower- or
#' upper-triangular; due to cache latency, this function will generally be
#' faster when `CF` is the lower Cholesky factor.
#' 
#' @param CF Cholesky factor (lower or upper) of dimension `p` \eqn{\times} `p`,
#' where e.g. `M = (CF)(CF')` is symmetric positive definite.
#' @param x A `p`-long vector such that `M - xx'` is positive definite
#' @param lower Boolean, dictating whether `CF` is the lower Cholesky factor
#' (`TRUE`, default) or upper Cholesky factor (`FASLE`). Also dictates whether
#' the returned Choleskyfactor is lower- or upper-triangular.
#'
#' @returns Cholesky factor of `M - xx'`
#' @export
#' @import Rcpp, RcppArmadillo
#' 
#' @examples
choldown = function(CF, x, lower = TRUE){
  
  p = length(x)
  R = CF
  R[1] = R[1] + 0
  x0 = x
  x0[1] = x0[1] + 0
  if(lower) choldownL_Rcpp(R, x0) else choldownU_Rcpp(R, x0)
  R
  
}