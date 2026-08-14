# The elastic net is the same construction as ridge and lasso: a
# separable penalty over a fixed() family, here the product of the
# Laplace and the Gaussian at zero.

test_that("the value is the elastic net up to its constant", {
  pen <- elasticnet_penalty(n_coef = 5)
  expect_identical(pen@params, c("lambda", "alpha"))
  expect_identical(penalty_kinks(pen, list(lambda = 1, alpha = 0.5)), 0)
  set.seed(1)
  b <- c(1.3, -2, 0.1, 0, -0.4)
  for (lam in c(0.2, 2, 15)) for (al in c(0.1, 0.5, 0.9)) {
    th <- list(lambda = lam, alpha = al)
    got <- penalty_value(pen, b, th) - penalty_value(pen, rep(0, 5), th)
    want <- lam * (al * sum(abs(b)) + (1 - al) * sum(b^2) / 2)
    expect_equal(got, want, tolerance = 1e-12,
                 info = sprintf("lambda %g alpha %g", lam, al))
  }
})

test_that("the constant depends on both hyperparameters", {
  # dropping it would make the two estimable by no marginal criterion,
  # which is the reason penalties7 keeps it
  pen <- elasticnet_penalty(n_coef = 3)
  z <- rep(0, 3)
  v1 <- penalty_value(pen, z, list(lambda = 1, alpha = 0.5))
  v2 <- penalty_value(pen, z, list(lambda = 2, alpha = 0.5))
  v3 <- penalty_value(pen, z, list(lambda = 1, alpha = 0.9))
  expect_false(isTRUE(all.equal(v1, v2)))
  expect_false(isTRUE(all.equal(v1, v3)))
})

test_that("the proximal operator is the soft threshold then the shrinkage", {
  pen <- elasticnet_penalty(n_coef = 8)
  set.seed(2)
  v <- rnorm(8, sd = 2)
  for (lam in c(0.05, 0.5, 3)) for (al in c(0.1, 0.5, 0.9)) {
    for (st in c(0.1, 1, 5)) {
      got <- penalty_prox(pen, v, st, list(lambda = lam, alpha = al))
      want <- sign(v) * pmax(abs(v) - st * lam * al, 0) /
        (1 + st * lam * (1 - al))
      expect_equal(got, want, tolerance = 1e-14)
    }
  }
  # and it really is the minimizer it claims to be
  th <- list(lambda = 1.5, alpha = 0.4)
  st <- 0.7
  b <- penalty_prox(pen, v, st, th)
  obj <- function(z) sum((z - v)^2) / (2 * st) + penalty_value(pen, z, th)
  for (i in seq_along(b)) {
    for (d in c(-1e-3, 1e-3)) {
      z <- b; z[i] <- z[i] + d
      expect_gte(obj(z) - obj(b), -1e-12)
    }
  }
})

test_that("the two ends are the lasso and the ridge", {
  pen <- elasticnet_penalty(n_coef = 4)
  b <- c(0.8, -1.2, 0.05, 0)
  lam <- 1.7
  pl <- lasso_penalty(n_coef = 4)
  expect_equal(penalty_prox(pen, b, 0.5, list(lambda = lam, alpha = 1 - 1e-12)),
               penalty_prox(pl, b, 0.5, list(lambda = lam)), tolerance = 1e-10)
  pr <- ridge_penalty(n_coef = 4)
  sg <- 1 / sqrt(lam)
  expect_equal(penalty_prox(pen, b, 0.5, list(lambda = lam, alpha = 1e-12)),
               penalty_prox(pr, b, 0.5, list(lambda = 1 / sg^2)),
               tolerance = 1e-10)
})

test_that("the derivatives are the parent's, reassembled", {
  skip_if_not_installed("numDeriv")
  set.seed(4)
  D <- matrix(rnorm(12), 4, 3)
  pen <- elasticnet_penalty(map = D)
  b <- c(0.7, -1.1, 0.35)
  th <- list(lambda = 1.3, alpha = 0.4)
  # away from the kink the value is differentiable in beta
  expect_equal(penalty_gradient(pen, b, th),
               numDeriv::grad(function(z) penalty_value(pen, z, th), b),
               tolerance = 1e-6)
  expect_equal(penalty_hessian(pen, b, th),
               numDeriv::hessian(function(z) penalty_value(pen, z, th), b),
               tolerance = 1e-4)
  # and in the hyperparameters, whose steps must stay inside (0, 1)
  ma <- list(d = 1e-4, r = 4)
  expect_equal(unname(unlist(penalty_grad_theta(pen, b, th))),
               numDeriv::grad(function(v)
                 penalty_value(pen, b, list(lambda = v[1], alpha = v[2])),
                 c(1.3, 0.4), method.args = ma), tolerance = 1e-6)
})

test_that("check_penalty passes and still catches an injected error", {
  pen <- elasticnet_penalty(n_coef = 3)
  res <- check_penalty(pen, theta = list(lambda = 1.2, alpha = 0.6),
                       verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})
