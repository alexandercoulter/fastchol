#' Cholesky Row/Column Addition
#'
#' @description
#' This function calculates the lower/upper Cholesky factor of positive definite
#' matrix `M`, provided the Cholesky factor of `M[-k, -k]` is known. That is,
#' this function inverts `choldrop`, re-inserting the `k`^th row/column from
#' `M`.
#' 
#' @param CF Cholesky factor, either lower-triangular (in which case set
#' `lower = TRUE`) or upper-triangular (in which case set `upper = FALSE`)
#' @param z The "missing" row/column of `M`, i.e. any `p + 1`-long vector which,
#' upon being inserted to be the `k`-th row/column in expanded matrix `M`,
#' retains positive definiteness `M > 0`
#' @param k The row/column index which `z`'s insertion will correspond to,
#' from `1` to `p + 1` inclusive
#' @param lower Boolean, dictating whether `CF` is the lower Cholesky factor
#' (`TRUE`, default) or upper Cholesky factor (`FASLE`). Also dictates whether
#' the returned Choleskyfactor is lower- or upper-triangular.
#'
#' @returns Cholesky factor of matrix `M`, where `M` is constructed from
#' `M[-k, -k]` and `z` inserted at the `k`-th row/column index.
#' @export
#' 
#' @examples
#' # Generate positive definite matrix
#' p = 20
#' M = diag(p) + tcrossprod(matrix(rnorm(p * 2), p, 2))
#' 
#' # Take out a row/column of M
#' k = 6
#' Mk = M[-k, -k]
#' x = M[ , k]
#' 
#' # Calculate Cholesky factor of Mk
#' Uk = chol(Mk)
#' 
#' # Calculate Cholesky factor adding back in the k-th row/column
#' L = choladd(t(Uk), x, k, lower = TRUE)
#' 
#' # Check against Cholesky factor of Mk
#' max(abs(L - t(chol(M))))
choladd = function(CF, z, k, lower = TRUE){
  
  p1 = length(z)
  R = matrix(0, p1, p1)
  if(lower) choladdL_Rcpp(R, CF, z, k) else choladdU_Rcpp(R, CF, z, k)
  R
  
}