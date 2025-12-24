#' Upper Cholesky Rank-1 Downdate
#'
#' @description
#' This function calculates the rank-1 downdate of upper Cholesky factor `U`,
#' i.e. the upper Cholesky factor of matrix `U'U - xx'` (provided that matrix is
#' also positive definite). Here, `U` is a square `p` \eqn{\times} `p`
#' upper-triangular matrix with positive diagonal entries, and `x` is any vector
#' of length `p` such that `U'U - xx' > 0`. This function applies hyperbolic
#' rotations to triangularize the joint matrix `(U x)`. Due to cache latency,
#' this function will generally be slower than `choldownL`.
#' 
#' @param U Upper Cholesky factor of dimension `p` \eqn{\times} `p`
#' @param x A `p`-long vector such that `U'U - xx' > 0`
#'
#' @returns Upper Cholesky factor of `U'U - xx'`
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
choldownU = function(U, x){
  
  p = length(x)
  U0 = U
  U0[1] = U0[1] + 0
  x0 = x
  x0[1] = x0[1] + 0
  choldownU_Rcpp(U0, x0)
  U0
  
}