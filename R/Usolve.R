#' Upper Triangular Solver
#'
#' @description
#' This function is a custom implementation of the default execution of 
#' `backsolve`, evaluating \eqn{U^{-1}x} where `U` is upper-triangular.
#' 
#' @param U Upper-triangular matrix of dimension `p` \eqn{\times} `p`
#' @param x Any vector of length `p`.
#'
#' @returns
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
Usolve = function(U, x){
  
  x0 = x
  x0[1] = x0[1] + 0
  Usolve_Rcpp(U, x0)
  x0
  
}