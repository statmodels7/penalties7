# The three branches, each against a route that shares no code with it.

test_that("the quadratic penalty passes its battery and knows its rank", {
  P <- crossprod(diff(diag(5), differences = 2))
  pen <- quadratic_penalty(P)
  res <- check_penalty(pen, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
  expect_identical(penalty_rank(pen), 3L)
  expect_false(is_proper(pen))
  expect_true(is_proper(quadratic_penalty(diag(4))))
  # the null space of a second-difference penalty is the linear polynomials
  nb <- penalty_null_basis(pen)
  expect_identical(ncol(nb), 2L)
  expect_lt(max(abs(crossprod(P, nb))), 1e-10)
})

test_that("the separable branch passes for ridge, lasso and the t prior", {
  for (pen in list(ridge_penalty(n_coef = 4), heavy_penalty(n_coef = 4))) {
    th <- if (identical(pen@params, "sigma")) list(sigma = 1.4) else
      list(sigma = 1.4, nu = 6)
    res <- check_penalty(pen, theta = th, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(pen@penalty_name,
                             paste(res$check[res$status != "OK"],
                                   collapse = ", ")))
  }
  # the lasso away from its declared kink
  pen <- lasso_penalty(n_coef = 3)
  expect_identical(penalty_kinks(pen, list(b = 1)), 0)
  res <- check_penalty(pen, beta = c(0.7, -1.2, 2.1), theta = list(b = 1),
                       verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("ridge exists twice and the two constructions agree exactly", {
  # quadratic_penalty(I) at lambda = 1/sigma^2 IS the gaussian at scale
  # sigma, up to the same constant: two implementations of one object, so
  # no tolerance has to be chosen beyond machine precision
  q <- 3
  beta <- c(0.4, -1.1, 2.2)
  sigma <- 1.7
  a <- quadratic_penalty(diag(q))
  b <- ridge_penalty(n_coef = q)
  tha <- list(lambda = 1 / sigma^2)
  thb <- list(sigma = sigma)
  expect_equal(penalty_value(a, beta, tha), penalty_value(b, beta, thb),
               tolerance = 1e-14)
  expect_equal(penalty_gradient(a, beta, tha), penalty_gradient(b, beta, thb),
               tolerance = 1e-14)
  expect_equal(penalty_hessian(a, beta, tha), penalty_hessian(b, beta, thb),
               tolerance = 1e-14)
})

test_that("SCAD and MCP match the literature and pass the battery", {
  pen <- scad_penalty(n_coef = 6)
  th <- list(lambda = 1, a = 3.7)
  # rho' typed independently from Fan-Li, on a grid inside the regions
  t <- c(0.3, 0.8, 1.5, 3.0, 4.5, 8.0)
  d1 <- penalty_gradient(pen, t, th)
  ref <- ifelse(t <= 1, 1, ifelse(t <= 3.7, (3.7 - t) / 2.7, 0))
  expect_equal(d1, ref, tolerance = 1e-14)
  expect_equal(penalty_value(pen, 0 * t, th), 0)
  # constant beyond the shoulder
  expect_equal(penalty_value(pen, 10, th), penalty_value(pen, 100, th))
  res <- check_penalty(pen, beta = t, theta = th, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))

  pen <- mcp_penalty(n_coef = 4)
  th <- list(lambda = 1, gamma = 3)
  t <- c(0.4, 1.2, 2.4, 5.0)
  d1 <- penalty_gradient(pen, t, th)
  expect_equal(d1, pmax(1 - t / 3, 0) * 1, tolerance = 1e-14)
  res <- check_penalty(pen, beta = t, theta = th, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("the map enters every derivative the same way", {
  D <- diff(diag(4))
  pen <- lasso_penalty(map = D)   # the fused lasso
  th <- list(b = 0.8)
  beta <- c(0.5, 1.4, -0.3, 2.0)
  res <- check_penalty(pen, beta = beta, theta = th, verbose = FALSE)
  expect_true(all(res$status == "OK"))
  # gradient through the map equals D' applied to the identity-map gradient
  flat <- lasso_penalty(n_coef = 3)
  g0 <- penalty_gradient(flat, as.numeric(D %*% beta), th)
  expect_equal(penalty_gradient(pen, beta, th),
               as.numeric(crossprod(D, g0)), tolerance = 1e-14)
})

test_that("the link scale is the chain rule on the parameter scale", {
  skip_if_not_installed("numDeriv")
  pen <- quadratic_penalty(diag(3))
  beta <- c(1, -0.5, 0.2)
  th <- list(lambda = 2)
  gl <- penalty_grad_theta(pen, beta, th, scale = "link")[[1]]
  ref <- numDeriv::grad(function(e) {
    penalty_value(pen, beta, list(lambda = exp(e)))
  }, log(2))
  expect_equal(gl, ref, tolerance = 1e-8)
  hl <- penalty_hess_theta(pen, beta, th, scale = "link")[[1]]
  refh <- numDeriv::hessian(function(e) {
    penalty_value(pen, beta, list(lambda = exp(e)))
  }, log(2))
  expect_equal(hl, as.numeric(refh), tolerance = 1e-6)
})

test_that("an injected defect is caught and a clean reference still passes", {
  pen <- quadratic_penalty(diag(3))
  res <- check_penalty(pen, verbose = FALSE)
  expect_true(all(res$status == "OK"))
  # a 5% wrong gradient must fail: probe through a corrupted subclass method
  Bad <- S7::new_class("Bad", parent = QuadraticPenalty)
  S7::method(penalty_gradient, Bad) <- function(pen, beta, theta, ...) {
    1.05 * theta[[1]] * as.numeric(pen@DPD %*% beta)
  }
  bad <- do.call(Bad, S7::props(quadratic_penalty(diag(3))))
  res <- check_penalty(bad, verbose = FALSE)
  expect_true(any(res$status == "FAILED"))
})

test_that("the base class refuses the marginal pieces off the quadratic branch", {
  pen <- ridge_penalty(n_coef = 2)
  expect_false(is_quadratic(pen))
  expect_error(penalty_matrix(pen, list(sigma = 1)), "quadratic")
  expect_error(penalty_rank(pen), "quadratic")
  expect_error(penalty_logpdet(pen, list(sigma = 1)), "quadratic")
})


test_that("the structured prior passes its battery on three structures", {
  skip_if_not_installed("numDeriv")
  for (s in list(parameters7::log_cholesky(3),
                 parameters7::ar1(4, role = "precision"),
                 parameters7::compound_symmetry(3, role = "precision"))) {
    pen <- structured_penalty(s)
    set.seed(5)
    th <- stats::setNames(as.list(round(rnorm(s@n_free, sd = 0.3), 2)),
                          s@free_names)
    res <- check_penalty(pen, theta = th, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(pen@penalty_name,
                             paste(res$check[res$status != "OK"],
                                   collapse = ", ")))
  }
})


test_that("the structured prior at the identity precision is the plain ridge", {
  # log_cholesky at a zero free vector is the identity matrix, so the prior
  # is quadratic_penalty(diag(q)) at lambda = 1: two constructions of one
  # object, compared with no tolerance to choose
  q <- 3
  beta <- c(0.4, -1.1, 2.2)
  s <- parameters7::log_cholesky(q)
  pen <- structured_penalty(s)
  th <- stats::setNames(as.list(rep(0, s@n_free)), s@free_names)
  ref <- quadratic_penalty(diag(q))
  expect_equal(penalty_value(pen, beta, th),
               penalty_value(ref, beta, list(lambda = 1)), tolerance = 1e-14)
  expect_equal(penalty_gradient(pen, beta, th),
               penalty_gradient(ref, beta, list(lambda = 1)),
               tolerance = 1e-14)
  expect_equal(unname(penalty_hessian(pen, beta, th)),
               penalty_hessian(ref, beta, list(lambda = 1)),
               tolerance = 1e-14)
  expect_true(is_quadratic(pen))
  expect_true(is_proper(pen))
  expect_identical(penalty_rank(pen), 3L)
  expect_equal(penalty_logpdet(pen, th)$value, 0, tolerance = 1e-14)
})
