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
#' 
#' @examples
#' # Generate positive definite matrix
#' p = 20
#' M = diag(p) + tcrossprod(matrix(rnorm(p * 2), p, 2))
#' 
#' # Calculate Cholesky factors
#' L = cholRcpp(M, lower = TRUE)
#' U = cholRcpp(M, lower = FALSE)
cholRcpp = function(M, lower = TRUE){
  
  p = ncol(M)
  R = matrix(0, p, p)
  cholRcpp_Rcpp(R, M, lower)
  R
  
}