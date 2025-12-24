#' Cholesky Row/Column Drop
#'
#' @description
#' This function calculates the Cholesky factor of matrix `M[-k, -k]`,
#' i.e. of positive definite `M` with the \eqn{k^{th}} row/column removed,
#' provided Cholesky factor `CF` (e.g. `M = (CF)(CF')`) is known. The user may
#' specify whether `CF` is lower- or upper-triangular; due to cache latency,
#' this function will generally be faster when `CF` is the lower Cholesky
#' factor.
#' 
#' @param CF Cholesky factor of dimension `p` \eqn{\times} `p`, such that e.g.
#' `(CF)(CF') = M`
#' @param k An index between `1` and `p` inclusive, corresponding to the
#' row/column of `M` to be dropped
#' @param lower Boolean, dictating whether `CF` is the lower Cholesky factor
#' (`TRUE`, default) or upper Cholesky factor (`FASLE`). Also dictates whether
#' the returned Cholesky factor is lower- or upper-triangular.
#'
#' @returns Cholesky factor of `M[-k, -k]`, i.e. of `M` with the
#' \eqn{k^{th}} row/column removed
#' @export
#' @import Rcpp, RcppArmadillo
#' 
#' @examples
choldrop = function(CF, k, lower = TRUE){
  
  p1 = ncol(CF) - 1
  R = matrix(0, p1, p1)
  if(lower) choldropL_Rcpp(R, CF, k) else choldropL_Rcpp(R, CF, k)
  R
  
}