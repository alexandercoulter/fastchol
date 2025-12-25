#' Triangular Solver
#'
#' @description
#' This function is a custom implementation of the default execution of 
#' `forwardsolve` (when `lower = TRUE`) or `backsolve` (when `lower = FALSE`),
#' evaluating \eqn{(CF)^{-1}x}.
#'
#' @param CF Triangular matrix of dimension `p` \eqn{\times} `p`
#' @param x Any vector of length `p`.
#' @param lower Boolean, dictating whether `CF` is lower-triangular
#' (`TRUE`, default) or upper-triangular (`FASLE`).
#'
#' @returns
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
trisolve = function(CF, x, lower = TRUE){
  
  x0 = x
  x0[1] = x0[1] + 0
  if(lower) Lsolve_Rcpp(CF, x0) else Usolve_Rcpp(CF, x0)
  x0
  
}