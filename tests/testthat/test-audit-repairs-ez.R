# The off-diagonal rows of check_penalty() and the validation of theta in the
# beta_quadratic() generic.

test_that("check_penalty() fails a wrong off-diagonal hyperparameter Hessian", {
  skip_if_not_installed("numDeriv")
  # SCAD with its mixed entry lambda_a shifted: the diagonal entries, the
  # gradients and the mixed blocks are those of the shipped class.
  BadMixed <- S7::new_class("BadMixedScad", parent = ScadPenalty)
  bad <- do.call(BadMixed, S7::props(scad_penalty(n_coef = 4)))
  S7::method(penalty_hess_theta, BadMixed) <- function(
      pen, beta, theta, scale = c("parameter", "link"), ...) {
    h <- S7::method(penalty_hess_theta, ScadPenalty)(pen, beta, theta, ...)
    h$lambda_a <- h$lambda_a + 0.01
    h
  }
  res <- check_penalty(bad, verbose = FALSE)
  row <- res[res$check == "hess_theta[lambda_a] vs numDeriv", ]
  expect_identical(nrow(row), 1L)
  expect_identical(row$status, "FAILED")
  expect_gt(row$max_error, 1e-3)
  expect_true(all(res$status[res$check != row$check] == "OK"))
})

test_that("check_penalty() has one passing row per pair of hyperparameters", {
  skip_if_not_installed("numDeriv")
  d1 <- crossprod(diff(diag(5)))
  d2 <- crossprod(diff(diag(5), differences = 2))
  # each penalty with its number of rows outside the off-diagonal ones
  cases <- list(
    scad = list(scad_penalty(n_coef = 4), 8L),
    mcp = list(mcp_penalty(n_coef = 4), 8L),
    elasticnet = list(elasticnet_penalty(n_coef = 4), 8L),
    heavy = list(heavy_penalty(n_coef = 4), 8L),
    additive2 = list(additive_penalty(list(d1, d2)), 13L),
    additive3 = list(additive_penalty(list(d1, diag(5), d2)), 15L),
    structured = list(structured_penalty(parameters7::ar1(4)), 12L),
    cholesky = list(structured_penalty(parameters7::log_cholesky(3)), 24L)
  )
  for (nm in names(cases)) {
    pen <- cases[[nm]][[1L]]
    res <- check_penalty(pen, verbose = FALSE)
    prs <- utils::combn(pen@params, 2L)
    keys <- paste0("hess_theta[", prs[1L, ], "_", prs[2L, ], "] vs numDeriv")
    expect_true(all(keys %in% res$check), label = nm)
    expect_identical(nrow(res), cases[[nm]][[2L]] + length(keys), label = nm)
    expect_true(all(res$status == "OK"), label = nm)
    expect_true(all(res$max_error[res$check %in% keys] < 1e-8), label = nm)
  }
})

test_that("beta_quadratic() checks theta on every branch", {
  q <- quadratic_penalty(diag(3))
  expect_error(beta_quadratic(q, list(lambda = -1)), "open interval")
  expect_error(beta_quadratic(q, list()), "lacks the parameter")
  a <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
  expect_error(beta_quadratic(a, list(lambda1 = 1)), "lacks the parameter")
  s <- structured_penalty(parameters7::log_cholesky(2))
  expect_error(beta_quadratic(s, list(log_L1 = 0, log_L2 = 0, L2.1 = Inf)),
               "open interval")
  k <- scad_penalty(n_coef = 3)
  expect_error(beta_quadratic(k, list(lambda = 1, a = 1.5)), "open interval")

  # a valid theta gives the same answers as before, as a list or a vector
  expect_true(beta_quadratic(q, list(lambda = 1)))
  expect_true(beta_quadratic(q, c(lambda = 1)))
  expect_true(beta_quadratic(a, list(lambda1 = 2, lambda2 = 0.5)))
  expect_true(beta_quadratic(s, list(log_L1 = 0, log_L2 = 0, L2.1 = 0.3)))
  expect_false(beta_quadratic(heavy_penalty(n_coef = 3),
                              list(sigma = 1, nu = 4)))
  expect_false(beta_quadratic(k, list(lambda = 1, a = 3.7)))
})
