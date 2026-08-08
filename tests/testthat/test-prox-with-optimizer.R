# The operator in the use it exists for: handed to a proximal method as the
# prox of the non-smooth part.
#
# The test lives here rather than in optimizers7 because the dependency runs
# this way: penalties7 sits above optimizers7 (through distributions7), so it
# may name it, while optimizers7 suggesting penalties7 would close a cycle --
# which is exactly what CI refused when the test was written on that side.

test_that("a penalty's operator drives a proximal run to the lasso solution", {
  skip_if_not_installed("optimizers7")

  set.seed(1)
  n <- 200
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  b0 <- c(2, -1.5, 0, 0, 0.8, 0)
  y <- as.numeric(X %*% b0 + rnorm(n))
  lambda <- 0.4

  fn <- function(b) sum((y - X %*% b)^2) / (2 * n)
  gr <- function(b) -as.numeric(crossprod(X, y - X %*% b)) / n

  pen <- lasso_penalty(n_coef = p)
  th <- list(lambda = lambda)
  fit <- optimizers7::minimize(
    optimizers7::prox_grad(
      prox = function(v, t) penalty_prox(pen, v, t, th),
      g = function(b) penalty_value(pen, b, th),
      criterion = optimizers7::crit_grad(1e-8)),
    fn = fn, gr = gr, par = rep(0, p))

  expect_true(fit@converged)

  # the stationarity conditions of the problem, written out here rather than
  # read from the run
  g <- gr(fit@par)
  on <- abs(fit@par) > 1e-8
  viol <- max(
    if (any(on)) max(abs(g[on] + lambda * sign(fit@par[on]))) else 0,
    if (any(!on)) max(pmax(abs(g[!on]) - lambda, 0)) else 0)
  expect_lt(viol, 1e-7)

  # the penalty keeps its normalizing constant, so the objective differs from
  # the bare lasso by that constant and the MINIMIZER does not
  bare <- optimizers7::minimize(
    optimizers7::prox_grad(
      prox = function(v, t) sign(v) * pmax(abs(v) - t * lambda, 0),
      g = function(b) lambda * sum(abs(b)),
      criterion = optimizers7::crit_grad(1e-8)),
    fn = fn, gr = gr, par = rep(0, p))
  expect_equal(fit@par, bare@par, tolerance = 1e-6)
})

test_that("the SCAD operator drives a run to a stationary point of its own objective", {
  skip_if_not_installed("optimizers7")

  set.seed(4)
  n <- 300
  p <- 8
  X <- matrix(rnorm(n * p), n, p)
  b0 <- c(3, -2, rep(0, p - 3), 1)
  y <- as.numeric(X %*% b0 + rnorm(n))

  fn <- function(b) sum((y - X %*% b)^2) / (2 * n)
  gr <- function(b) -as.numeric(crossprod(X, y - X %*% b)) / n
  pen <- scad_penalty(n_coef = p)
  th <- list(lambda = 0.15, a = 3.7)

  fit <- optimizers7::minimize(
    optimizers7::prox_grad(
      prox = function(v, t) penalty_prox(pen, v, t, th),
      g = function(b) penalty_value(pen, b, th),
      step = 0.5, criterion = optimizers7::crit_grad(1e-8)),
    fn = fn, gr = gr, par = rep(0, p))
  expect_true(fit@converged)

  # away from the kinks SCAD is differentiable, so the total gradient must
  # vanish there; the large coefficients sit past a*lambda, where the penalty
  # is flat and the fit is the unpenalized one
  b <- fit@par
  free <- abs(b) > th$a * th$lambda
  expect_true(any(free))
  expect_lt(max(abs((gr(b) + penalty_gradient(pen, b, th))[free])), 1e-6)
})
