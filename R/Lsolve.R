#' Lower Triangular Solver
#'
#' @description
#' This function is a custom implementation of the default execution of 
#' `forwardsolve`, evaluating \eqn{L^{-1}x} where `L` is lower-triangular.
#'
#' @param L Lower-triangular matrix of dimension `p` \eqn{\times} `p`
#' @param x Any vector of length `p`.
#'
#' @returns
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
Lsolve = function(L, x){
  
  x0 = x
  x0[1] = x0[1] + 0
  Lsolve_Rcpp(L, x0)
  x0
  
}