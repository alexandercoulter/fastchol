#' Cholesky Decomposition with Rcpp
#'
#' @description
#' This function calculates the Cholesky decomposition of matrix `M` with the
#' Rcpp library implementation, which permits calculating either the lower or
#' upper Cholesky factors. Calculating the lower Cholesky factor tends to be
#' twice as fast, so it is the default.
#' 
#' @param M Any square positive definite matrix
#' @param lower Boolean, dictating whether the lower Cholesky factor (`TRUE`,
#' default) or upper Cholesky factor (`FASLE`) should be returned
#'
#' @returns The lower (or upper) Cholesky factor of `M`
#' @export
#' @import Rcpp, RcppArmadillo
#' 
#' @examples
cholRcpp = function(M, lower = TRUE){
  
  p = ncol(M)
  R = matrix(0, p, p)
  cholRcpp_Rcpp(R, M, lower)
  R
  
}