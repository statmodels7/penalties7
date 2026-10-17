# The length of the coefficients, the shape of the hyperparameters and the
# validation inside beta_quadratic().

# One penalty per branch, with a valid theta and its number of coefficients.
g1_branches <- function() {
  list(
    quadratic = list(pen = quadratic_penalty(diag(3)),
                     theta = list(lambda = 2), q = 3L),
    mapped = list(pen = quadratic_penalty(diag(2), map = diff(diag(3))),
                  theta = list(lambda = 2), q = 3L),
    lasso = list(pen = lasso_penalty(n_coef = 3),
                 theta = list(lambda = 1), q = 3L),
    heavy = list(pen = heavy_penalty(n_coef = 3),
                 theta = list(sigma = 1, nu = 4), q = 3L),
    scad = list(pen = scad_penalty(n_coef = 3),
                theta = list(lambda = 1, a = 3.7), q = 3L),
    mcp = list(pen = mcp_penalty(n_coef = 3),
               theta = list(lambda = 1, gamma = 3), q = 3L),
    additive = list(pen = additive_penalty(list(diag(3), diag(c(1, 1, 0)))),
                    theta = list(lambda1 = 2, lambda2 = 0.5), q = 3L),
    structured = list(pen = structured_penalty(parameters7::log_cholesky(2)),
                      theta = list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3),
                      q = 2L)
  )
}

# Every exported generic that receives the coefficients, called with a
# coefficient vector b and directions of the right length.
g1_calls <- function(pen, b, theta, q) {
  v <- seq_len(q) / 10
  w <- rev(v)
  list(
    penalty_value = function() penalty_value(pen, b, theta),
    penalty_gradient = function() penalty_gradient(pen, b, theta),
    penalty_hessian = function() penalty_hessian(pen, b, theta),
    penalty_grad_theta = function() penalty_grad_theta(pen, b, theta),
    penalty_hess_theta = function() penalty_hess_theta(pen, b, theta),
    penalty_cross = function() penalty_cross(pen, b, theta),
    penalty_dhessian = function() penalty_dhessian(pen, b, theta),
    penalty_d2hessian = function() penalty_d2hessian(pen, b, theta),
    penalty_dcross = function() penalty_dcross(pen, b, theta),
    penalty_dhessian_beta = function() penalty_dhessian_beta(pen, b, theta, v),
    penalty_d2hessian_beta = function()
      penalty_d2hessian_beta(pen, b, theta, v, w),
    penalty_dhessian_beta_theta = function()
      penalty_dhessian_beta_theta(pen, b, theta, v)
  )
}

test_that("every generic that receives beta checks its length", {
  for (nm in names(g1_branches())) {
    br <- g1_branches()[[nm]]
    msg <- sprintf("'beta' must have length %d", br$q)
    for (b in list(seq_len(br$q + 2L) / 3, seq_len(br$q - 1L) / 3)) {
      calls <- g1_calls(br$pen, b, br$theta, br$q)
      for (g in names(calls)) {
        expect_error(calls[[g]](), msg, fixed = TRUE,
                     info = paste(nm, g, length(b)))
      }
    }
  }
  # the right length still works; the marginal generics exist only for a
  # branch without a kink
  smooth <- c("quadratic", "mapped", "heavy", "additive", "structured")
  first <- c("penalty_value", "penalty_gradient", "penalty_hessian",
             "penalty_grad_theta", "penalty_hess_theta", "penalty_cross")
  for (nm in names(g1_branches())) {
    br <- g1_branches()[[nm]]
    b <- c(0.7, -0.4, 1.3)[seq_len(br$q)]
    calls <- g1_calls(br$pen, b, br$theta, br$q)
    use <- if (nm %in% smooth) names(calls) else first
    for (g in use) expect_no_error(calls[[g]](), message = paste(nm, g))
    expect_length(penalty_gradient(br$pen, b, br$theta), br$q)
  }
  # an integer vector and a one-column matrix are still coerced
  pen <- lasso_penalty(n_coef = 3)
  b <- c(2L, -1L, 3L)
  ref <- -sum(log(dexp(abs(b), rate = 1.5) / 2))
  expect_equal(penalty_value(pen, b, list(lambda = 1.5)), ref,
               tolerance = 1e-14)
  expect_equal(penalty_value(pen, matrix(b, ncol = 1), list(lambda = 1.5)),
               ref, tolerance = 1e-14)
})

test_that("a hyperparameter that is not a single number is an error", {
  pen <- quadratic_penalty(diag(3))
  b <- c(1, -0.5, 0.25)
  for (lam in list(c(1, 2), numeric(0), NULL)) {
    th <- list(lambda = lam)
    expect_error(penalty_value(pen, b, th),
                 "Parameter 'lambda' must be a single number.", fixed = TRUE)
    expect_error(penalty_gradient(pen, b, th),
                 "Parameter 'lambda' must be a single number.", fixed = TRUE)
    expect_error(penalty_matrix(pen, th),
                 "Parameter 'lambda' must be a single number.", fixed = TRUE)
  }
  sep <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 2)
  expect_error(penalty_value(sep, c(1, -1), list(sigma = c(1, 2))),
               "Parameter 'sigma' must be a single number.", fixed = TRUE)
  expect_error(penalty_grad_theta(sep, c(1, -1), list(sigma = c(1, 2))),
               "Parameter 'sigma' must be a single number.", fixed = TRUE)
  expect_error(penalty_kinks(scad_penalty(), list(lambda = c(1, 2), a = 3.7)),
               "Parameter 'lambda' must be a single number.", fixed = TRUE)
  # the open-interval message is unchanged
  expect_error(penalty_value(pen, b, list(lambda = -1)),
               "open interval", fixed = TRUE)
  expect_error(penalty_value(pen, b, list(lambda = NA_real_)),
               "open interval", fixed = TRUE)
  # and a scalar works, in a list or a named vector, against the Gaussian
  # log-density written independently
  ref <- -sum(stats::dnorm(b, sd = 1 / sqrt(2), log = TRUE))
  expect_equal(penalty_value(pen, b, list(lambda = 2)), ref,
               tolerance = 1e-14)
  expect_equal(penalty_value(pen, b, c(lambda = 2)), ref, tolerance = 1e-14)
  expect_equal(penalty_value(sep, c(1, -1), list(sigma = 1.5)),
               -sum(stats::dnorm(c(1, -1), sd = 1.5, log = TRUE)),
               tolerance = 1e-14)
})

test_that("beta_quadratic() signals an invalid theta", {
  h <- heavy_penalty(n_coef = 2)
  expect_error(beta_quadratic(h, list(sigma = -1, nu = 4)),
               "Parameter 'sigma' must lie in the open interval", fixed = TRUE)
  expect_error(beta_quadratic(h, list(sigma = 1)),
               "'theta' lacks the parameter(s) nu;", fixed = TRUE)
  expect_error(beta_quadratic(h, list(sigma = c(1, 2), nu = 4)),
               "Parameter 'sigma' must be a single number.", fixed = TRUE)
  # a valid theta still gets the answer: a Student t prior is not quadratic in
  # the coefficients and a Gaussian one is
  expect_false(beta_quadratic(h, list(sigma = 1, nu = 4)))
  g <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 2)
  expect_true(beta_quadratic(g, list(sigma = 1)))
})

test_that("beta_quadratic() is TRUE for a parent answering one curvature", {
  # a Gaussian whose response curvature is returned as one number for all
  # the points, as a family with a constant curvature may: the log-density is
  # then quadratic in the response
  OneHess <- S7::new_class(
    "OneHess", parent = S7::S7_class(distributions7::gaussian1_distrib()))
  hy <- distributions7::distrib_hess_y
  S7::method(hy, OneHess) <- function(distrib, y, theta, ...) {
    -1 / theta[["sigma"]]^2
  }
  fake <- S7::convert(distributions7::gaussian1_distrib(), OneHess)
  expect_length(distributions7::distrib_hess_y(fake, c(-1, 0, 2),
                                               list(mu = 0, sigma = 2)), 1L)
  pen <- distrib_penalty(distributions7::fixed(fake, mu = 0), n_coef = 3)
  expect_true(beta_quadratic(pen, list(sigma = 2)))
  # and the movement of its Hessian with the coefficients is then zero
  expect_identical(penalty_dhessian_beta(pen, c(1, -1, 0.5), list(sigma = 2),
                                         c(0.1, 0.2, 0.3)),
                   matrix(0, 3, 3))
})
