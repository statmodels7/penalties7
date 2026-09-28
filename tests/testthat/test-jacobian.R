# A penalty written on D beta is a prior on beta only once the Jacobian of
# beta -> D beta is in its constant. These tests pin that the value is the
# negative log-density of beta under a diagonal map, and that it does not
# move when a coordinate is rescaled together with its map entry.

test_that("a quadratic penalty under a diagonal map carries the Jacobian", {
  d <- c(0.5, 2, 3.7, 0.1)
  b <- c(0.3, -1.2, 0.8, 4)
  P <- diag(4)
  mapped <- quadratic_penalty(P, map = Matrix::Diagonal(x = d))
  plain <- quadratic_penalty(P)
  for (lam in c(0.2, 1, 7)) {
    th <- list(lambda = lam)
    # the density of beta is the density of u = D beta times |det D|
    expect_equal(penalty_value(mapped, b, th),
                 penalty_value(plain, d * b, th) - sum(log(d)),
                 tolerance = 1e-12)
    # and it is the Gaussian density with precision lambda D'PD
    expect_equal(penalty_value(mapped, b, th),
                 -sum(stats::dnorm(b, sd = 1 / (sqrt(lam) * d), log = TRUE)),
                 tolerance = 1e-12)
  }
  expect_equal(mapped@logpdet_DPD, 2 * sum(log(d)), tolerance = 1e-12)
  expect_identical(penalty_rank(mapped), 4L)
})

test_that("a deficient quadratic penalty under a map reads pdet(D'PD)", {
  d <- c(1.5, 0.4, 2, 0.7, 3)
  P <- crossprod(diff(diag(5), differences = 2))
  mapped <- quadratic_penalty(P, map = Matrix::Diagonal(x = d))
  DPD <- diag(d) %*% P %*% diag(d)
  ev <- eigen(DPD, symmetric = TRUE)$values
  ev <- ev[ev > 1e-10 * max(ev)]
  expect_identical(penalty_rank(mapped), 3L)
  expect_equal(mapped@logpdet_DPD, sum(log(ev)), tolerance = 1e-10)
  lp <- penalty_logpdet(mapped, list(lambda = 2))
  expect_equal(lp$value, 3 * log(2) + sum(log(ev)), tolerance = 1e-10)
})

test_that("a separable penalty under a diagonal map carries the Jacobian", {
  d <- c(0.5, 2, 3.7)
  b <- c(0.3, -1.2, 0.8)
  mapped <- lasso_penalty(map = Matrix::Diagonal(x = d))
  plain <- lasso_penalty(n_coef = 3)
  th <- list(lambda = 1.3)
  expect_equal(penalty_value(mapped, b, th),
               penalty_value(plain, d * b, th) - sum(log(d)),
               tolerance = 1e-12)
  # no derivative moves: the term is constant in beta and in theta
  expect_equal(penalty_gradient(mapped, b, th),
               d * penalty_gradient(plain, d * b, th), tolerance = 1e-12)
  expect_equal(penalty_grad_theta(mapped, b, th)$lambda,
               penalty_grad_theta(plain, d * b, th)$lambda, tolerance = 1e-12)
})

test_that("the value does not move when a coordinate changes its units", {
  # A column measured in other units is x c, its standard deviation sd c, its
  # coefficient beta / c and the map entry d / c: D beta is unchanged, and so
  # must be the density of the data, hence the prior times the change of
  # variable. Here that reads: the value at (beta / c, d c) equals the value
  # at (beta, d) plus log(c), the Jacobian of the relabelling itself.
  d <- c(0.5, 2, 3.7)
  b <- c(0.3, -1.2, 0.8)
  cc <- c(10, 1, 0.01)
  for (mk in list(function(m) ridge_penalty(map = m),
                  function(m) lasso_penalty(map = m))) {
    a <- mk(Matrix::Diagonal(x = d))
    z <- mk(Matrix::Diagonal(x = d * cc))
    th <- list(lambda = 0.8)
    expect_equal(penalty_value(z, b / cc, th),
                 penalty_value(a, b, th) - sum(log(cc)), tolerance = 1e-12)
  }
})

test_that("a penalty without a map is unchanged", {
  P <- crossprod(diff(diag(6), differences = 2))
  pen <- quadratic_penalty(P)
  ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  keep <- ev > 1e-10 * max(ev)
  expect_identical(pen@logpdet_DPD, sum(log(ev[keep])))
  expect_identical(penalty_rank(pen), sum(keep))
  b <- c(1, 0, -1)
  expect_identical(penalty_value(lasso_penalty(n_coef = 3), b,
                                 list(lambda = 0.5)),
                   0.5 * sum(abs(b)) - 3 * log(0.5 / 2))
})
