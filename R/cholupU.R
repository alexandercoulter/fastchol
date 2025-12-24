#' Upper Cholesky Rank-1 Update
#'
#' @description
#' This function calculates the rank-1 update of upper Cholesky factor `U`, i.e.
#' the upper Cholesky factor of matrix `U'U + xx'`. Here, `U` is a square `p`
#' \eqn{\times} `p` upper-triangular matrix with positive diagonal entries, and
#' `x` is any arbitrary vector of length `p`. This function applies Givens
#' rotations to triangularize the joint matrix `(U x)`. Due to cache latency,
#' this function will generally be slower than `cholupL`.
#' 
#' @param U Upper Cholesky factor of dimension `p` \eqn{\times} `p`
#' @param x Any `p`-long vector
#'
#' @returns Upper Cholesky factor of `U'U + xx'`
#' @export
#' @import Rcpp, RcppArmadillo
#' 
#' @examples
cholupU = function(U, x){
  
  p = length(x)
  U0 = U
  U0[1] = U0[1] + 0
  x0 = x
  x0[1] = x0[1] + 0
  cholupU_Rcpp(U0, x0)
  U0
  
}