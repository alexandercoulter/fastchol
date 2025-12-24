#' Lower Cholesky Row/Column Addition
#'
#' @description
#' This function calculates the lower Cholesky factor of positive definite
#' matrix `M`, provided lower Cholesky factor of `M[-k, -k]` is known. That is,
#' this function inverts `choldropL`, re-inserting the `k`^th row/column from
#' `M`.
#' 
#' @param L Lower Cholesky factor of dimension `p` \eqn{\times} `p`, such that
#' `LL' = M[-k, -k]`
#' @param z The "missing" row/column of `M`, i.e. any `p + 1`-long vector which,
#' upon being inserted to be the `k`-th row/column in expanded matrix `M`,
#' retains positive definiteness `M > 0`
#' @param k The row/column index which `z`'s insertion will correspond to,
#' from `1` to `p + 1` inclusive
#'
#' @returns Lower Cholesky factor of matrix `M`, where `M` is constructed from
#' `M[-k, -k]` and `z` inserted at the `k`-th row/column index.
#' @export
#' @import Rcpp, RcppArmadillo
#' 
#' @examples
choladdL = function(L, z, k){
  
  p1 = length(z)
  L1 = matrix(0, p1, p1)
  choladdL_Rcpp(L1, L, z, k)
  L1
  
}