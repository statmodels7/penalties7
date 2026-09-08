# Drawing from a penalty read as a prior. Every check is against a route that
# shares no arithmetic with the draw: the parent's own distribution function,
# or the inverse of the precision matrix the Gaussian branches are written
# from. A draw that ignored its hyperparameter, or inverted a map the wrong
# way, would pass neither.

test_that("a separable draw follows the parent's own distribution", {
  set.seed(1)
  g <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 200)
  x <- as.numeric(replicate(40, penalty_draw(g, list(sigma = 0.4))))
  # the parent's cdf, which the draw does not read
  p <- distributions7::distrib_cdf(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    x, list(sigma = 0.4))
  expect_gt(stats::ks.test(p, "punif")$p.value, 1e-4)
  expect_equal(stats::sd(x), 0.4, tolerance = 0.02)

  # and the hyperparameter is read: the same call at another scale
  y <- as.numeric(replicate(40, penalty_draw(g, list(sigma = 1.6))))
  expect_equal(stats::sd(y), 1.6, tolerance = 0.02)
})

test_that("a Laplace prior gives Laplace effects, not Gaussian ones", {
  set.seed(2)
  l <- lasso_penalty(n_coef = 400)
  x <- as.numeric(replicate(20, penalty_draw(l, list(lambda = 2))))
  p <- distributions7::distrib_cdf(l@parent, x, list(lambda = 2))
  expect_gt(stats::ks.test(p, "punif")$p.value, 1e-4)
  # the negative control: it is NOT the Gaussian of the same variance
  q <- stats::pnorm(x, 0, stats::sd(x))
  expect_lt(stats::ks.test(q, "punif")$p.value, 1e-6)
})

test_that("a heavy-tailed prior is drawn at its own degrees of freedom", {
  set.seed(3)
  h <- heavy_penalty(n_coef = 300)
  x <- as.numeric(replicate(20, penalty_draw(h, list(sigma = 1, nu = 3))))
  p <- distributions7::distrib_cdf(h@parent, x, list(sigma = 1, nu = 3))
  expect_gt(stats::ks.test(p, "punif")$p.value, 1e-4)
})

test_that("a Gaussian branch draws at the covariance its precision implies", {
  set.seed(4)
  # a full-rank, non-diagonal precision, so the Cholesky route runs
  A <- matrix(c(2, 0.8, 0.3, 0.8, 1.5, -0.4, 0.3, -0.4, 1.2), 3, 3)
  pen <- quadratic_penalty(A)
  expect_true(is_proper(pen))
  x <- t(replicate(60000, penalty_draw(pen, list(lambda = 1.5))))
  expect_equal(stats::cov(x), solve(as.matrix(penalty_matrix(pen,
                                                             list(lambda = 1.5)))),
               tolerance = 0.03)
})

test_that("a diagonal precision draws coordinate by coordinate", {
  set.seed(5)
  pen <- quadratic_penalty(diag(c(1, 4, 9)))
  x <- t(replicate(20000, penalty_draw(pen, list(lambda = 2))))
  expect_equal(apply(x, 2, stats::sd), 1 / sqrt(2 * c(1, 4, 9)),
               tolerance = 0.03)
  # off-diagonal covariance is zero, so the draw is not correlated
  expect_lt(max(abs(stats::cov(x)[upper.tri(diag(3))])), 0.01)
})

test_that("a direction no prior covers comes back NA", {
  set.seed(6)
  # the Demmler-Reinsch shape: the linear column is unpenalized
  pen <- quadratic_penalty(diag(c(0, 1, 1, 1)))
  d <- penalty_draw(pen, list(lambda = 1))
  expect_true(is.na(d[[1L]]))
  expect_true(all(is.finite(d[-1L])))

  # SCAD and MCP are densities of nothing
  expect_true(all(is.na(penalty_draw(scad_penalty(n_coef = 4),
                                     list(lambda = 1, a = 3.7)))))
  expect_true(all(is.na(penalty_draw(mcp_penalty(n_coef = 4),
                                     list(lambda = 1, gamma = 3)))))

  # and a deficient non-diagonal precision, which is what a tensor smooth has
  P <- crossprod(diff(diag(5), differences = 2))
  expect_true(all(is.na(penalty_draw(quadratic_penalty(P),
                                     list(lambda = 1)))))
})

test_that("a diagonal map is inverted and a general one is refused", {
  set.seed(7)
  d <- c(0.5, 2, 4)
  m <- Matrix::Diagonal(x = d)
  g <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 3, map = m)
  x <- t(replicate(20000, penalty_draw(g, list(sigma = 1))))
  # the prior is on d * beta, so beta has standard deviation sigma / d
  expect_equal(apply(x, 2, stats::sd), 1 / d, tolerance = 0.03)

  # a general map leaves the coefficients underdetermined
  gen <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 3, map = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2))
  expect_true(all(is.na(penalty_draw(gen, list(sigma = 1)))))
})

test_that("a multivariate prior draws whole blocks", {
  set.seed(8)
  mv <- distrib_penalty(
    distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                          mu1 = 0, mu2 = 0), n_coef = 2 * 4000)
  th <- list(sigma_log_L1 = log(0.5), sigma_log_L2 = log(1.2),
             sigma_L2.1 = 0.7)
  x <- matrix(penalty_draw(mv, th), ncol = 2, byrow = TRUE)
  # against the covariance the same hyperparameters describe, built by the
  # structure rather than by the draw
  S <- distributions7::mv_sigma(mv@parent, th)
  expect_equal(stats::cov(x), as.matrix(S), tolerance = 0.06,
               ignore_attr = TRUE)
})

test_that("the draw consumes the caller's stream and nothing else", {
  pen <- ridge_penalty(n_coef = 5)
  set.seed(9)
  a <- penalty_draw(pen, list(lambda = 1))
  set.seed(9)
  b <- penalty_draw(pen, list(lambda = 1))
  expect_identical(a, b)
  # two consecutive draws differ, so no seed is set inside
  set.seed(9)
  expect_false(isTRUE(all.equal(penalty_draw(pen, list(lambda = 1)),
                                penalty_draw(pen, list(lambda = 1)))))
})
