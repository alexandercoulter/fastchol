#' Lower Cholesky Rank-1 Downdate
#'
#' @description
#' This function calculates the rank-1 downdate of lower Cholesky factor `L`,
#' i.e. the lower Cholesky factor of matrix `LL' - xx'` (provided that matrix is
#' also positive definite). Here, `L` is a square `p` \eqn{\times} `p`
#' lower-triangular matrix with positive diagonal entries, and `x` is any vector
#' of length `p` such that `LL' - xx' > 0`. This function applies hyperbolic
#' rotations to triangularize the joint matrix `(L x)`. Due to cache latency,
#' this function will generally be faster than `choldownU`.
#' 
#' @param L Lower Cholesky factor of dimension `p` \eqn{\times} `p`
#' @param x A `p`-long vector such that `LL' - xx' > 0`
#'
#' @returns Lower Cholesky factor of `LL' - xx'`
#' @export
#' @import Rcpp, RcppArmadillo
#' 
#' @examples
choldownL = function(L, x){
  
  p = length(x)
  L0 = L
  L0[1] = L0[1] + 0
  x0 = x
  x0[1] = x0[1] + 0
  choldownL_Rcpp(L0, x0)
  L0
  
}