#' Lower Cholesky Row/Column Drop
#'
#' @description
#' This function calculates the lower Cholesky factor of matrix `M[-k, -k]`,
#' i.e. of positive definite `M` with the \eqn{k^{th}} row/column removed,
#' provided lower Cholesky factor `L` (`M = LL'`) is known. Due to cache
#' latency, this function will generally be faster than `choldropU`.
#' 
#' @param L Lower Cholesky factor of dimension `p` \eqn{\times} `p`, such that
#' `LL' = M`
#' @param k An index between `1` and `p` inclusive, corresponding to the
#' row/column to be dropped
#'
#' @returns Lower Cholesky factor of `M[-k, -k]`, i.e. of `M` with the
#' \eqn{k^{th}} row/column removed
#' @export
#' @import Rcpp, RcppArmadillo
#' 
#' @examples
choldropL = function(L, k){
  
  p1 = ncol(L) - 1
  L0 = matrix(0, p1, p1)
  choldropL_Rcpp(L0, L, k)
  L0
  
}