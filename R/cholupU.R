#' Cholesky Rank-1 Update Using Upper Cholesky Factor U
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