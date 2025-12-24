#' Cholesky Rank-1 Update Using Lower Cholesky Factor L
#'
#' @param L Lower Cholesky factor of dimension `p` \eqn{\times} `p`
#' @param x Any `p`-long vector
#'
#' @returns Lower Cholesky factor of `LL' + xx'`
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
cholupL = function(L, x){
  
  p = length(x)
  L0 = L
  L0[1] = L0[1] + 0
  x0 = x
  x0[1] = x0[1] + 0
  cholupL_Rcpp(L0, x0)
  L0
  
}