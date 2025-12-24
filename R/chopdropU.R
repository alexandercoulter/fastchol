#' Uper Cholesky Row/Column Drop
#'
#' @description
#' This function calculates the upper Cholesky factor of matrix `M[-k, -k]`,
#' i.e. of positive definite `M` with the \eqn{k^{th}} row/column removed,
#' provided upper Cholesky factor `U` (`M = U'U`) is known. Due to cache
#' latency, this function will generally be slower than `choldropL`.
#' 
#' @param U Upper Cholesky factor of dimension `p` \eqn{\times} `p`, such that
#' `U'U = M`
#' @param k An index between `1` and `p` inclusive, corresponding to the
#' row/column to be dropped
#'
#' @returns Upper Cholesky factor of `M[-k, -k]`, i.e. of `M` with the
#' \eqn{k^{th}} row/column removed
#' @export
#' @import Rcpp, RcppArmadillo
#' 
#' @examples
choldropU = function(U, k){
  
  p1 = ncol(U) - 1
  U0 = matrix(0, p1, p1)
  choldropU_Rcpp(U0, U, k)
  R0
  
}