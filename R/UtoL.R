#' Fast U to L Transpose
#'
#' @description
#' This function copies the upper-triangular part of square matrix `U`
#' (including the main diagonal) into the lower-triangular portion of output
#' matrix `L`, other values zero. For `U` that is only upper-triangular, this is
#' effectively a faster transpose operation.
#' 
#' @param U Any square matrix such that its upper-triangular part (including
#' main diagonal) will be transposed to create a lower-triangular matrix.
#'
#' @returns A lower-triangular square matrix whose non-zero entries are the
#' upper-triangular part of `U`.
#' @export
#' @import Rcpp, RcppArmadillo
#'
#' @examples
UtoL = function(U){
  
  L = matrix(0, nrow(U), ncol(U))
  UtoL_Rcpp(L, U)
  L
  
}