test_that("a multivariate t prior's Hessian moves with the coefficients exactly", {
  f <- distributions7::fixed(distributions7::mvstudent_t1_distrib(2),
                             mu1 = 0, mu2 = 0)
  pen <- distrib_penalty(f, n_coef = 8)
  th <- list(sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.3,
             nu = 3)
  expect_false(beta_quadratic(pen, th))
  set.seed(1)
  b <- rnorm(8)
  v <- rnorm(8)
  w <- rnorm(8)
  h <- 1e-5

  D <- penalty_dhessian_beta(pen, b, th, v)
  num <- (penalty_hessian(pen, b + h * v, th) -
            penalty_hessian(pen, b - h * v, th)) / (2 * h)
  expect_lt(max(abs(D - num)), 1e-7)
  # block diagonal: the blocks are independent under the prior
  expect_identical(D[1:2, 3:8], matrix(0, 2, 6))

  D2 <- penalty_d2hessian_beta(pen, b, th, v, w)
  num <- (penalty_dhessian_beta(pen, b + h * w, th, v) -
            penalty_dhessian_beta(pen, b - h * w, th, v)) / (2 * h)
  expect_lt(max(abs(D2 - num)), 1e-7)

  Dt <- penalty_dhessian_beta_theta(pen, b, th, v)
  expect_named(Dt, names(th))
  for (m in names(th)) {
    tp <- th
    tp[[m]] <- tp[[m]] + h
    tm <- th
    tm[[m]] <- tm[[m]] - h
    num <- (penalty_dhessian_beta(pen, b, tp, v) -
              penalty_dhessian_beta(pen, b, tm, v)) / (2 * h)
    expect_lt(max(abs(Dt[[m]] - num)), 1e-7)
  }
})

test_that("a multivariate gaussian prior still answers zero without its parent", {
  g <- distrib_penalty(
    distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                          mu1 = 0, mu2 = 0), n_coef = 6)
  th <- list(sigma_log_L1 = 0.1, sigma_log_L2 = -0.2, sigma_L2.1 = 0.3)
  b <- c(1, -0.5, 0.3, 2, 0.1, -1)
  expect_identical(penalty_dhessian_beta(g, b, th, rev(b)), matrix(0, 6, 6))
})

test_that("dp_contract() contracts the trailing index against each direction", {
  set.seed(3)
  A <- array(rnorm(2 * 2 * 2 * 3), c(2, 2, 2, 3))
  v <- matrix(rnorm(6), 3, 2)
  out <- dp_contract(A, list(v))
  for (i in 1:3) {
    ref <- A[, , 1, i] * v[i, 1] + A[, , 2, i] * v[i, 2]
    expect_equal(out[, , i], ref, tolerance = 1e-14)
  }
})
