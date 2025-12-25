test_that("cholup accurately gives the Cholesky factor of rank-1 update", {
  # Set parameters
  p = 10
  
  # Generate positive definite matrix
  M = diag(p) + tcrossprod(matrix(rnorm(2 * p), p, 2))
  
  # Generate vector to update with
  x = rnorm(p)
  
  # Calculate rank-1 updated matrix
  Mx = M + tcrossprod(x)
  
  # Calculate Cholesky factors (upper, lower) of M, Mx
  UM = chol(M)
  LM = t(UM)
  UMx = chol(Mx)
  LMx = t(UMx)
  
  # Calculate rank-1 updates using cholup
  UMx_cholup = cholup(UM, x, lower = FALSE)
  LMx_cholup = cholup(LM, x, lower = TRUE)
  
  # Calculate error from Cholesky factors using chol(Mx)
  errors = abs(UMx - UMx_cholup) + abs(LMx - LMx_cholup)
  
  expect_true(
    max(errors) <= 1e-10
  )
})
