# What a marginal criterion asks of a penalty: the theta-derivatives of the
# coefficient Hessian and of the mixed block.

# numDeriv differentiates penalty_hessian and penalty_cross themselves, so the
# reference shares no arithmetic with any of the methods under test.
num_dhessian <- function(pen, beta, theta) {
  p <- pen@params
  v0 <- vapply(p, function(q) as.numeric(theta[[q]]), numeric(1))
  lapply(stats::setNames(seq_along(p), p), function(m) {
    f <- function(x) {
      th <- as.list(v0)
      th[[m]] <- x
      as.numeric(penalty_hessian(pen, beta, th))
    }
    k <- length(beta)
    matrix(numDeriv::jacobian(f, v0[[m]]), k, k)
  })
}

num_dcross <- function(pen, beta, theta) {
  p <- pen@params
  v0 <- vapply(p, function(q) as.numeric(theta[[q]]), numeric(1))
  prs <- penalties7:::ptheta_pairs(p)
  stats::setNames(lapply(names(prs), function(nm) {
    pr <- prs[[nm]]
    f <- function(x) {
      th <- as.list(v0)
      th[[pr[1]]] <- x
      as.numeric(penalty_cross(pen, beta, th)[[pr[2]]])
    }
    as.numeric(numDeriv::jacobian(f, v0[[pr[1]]]))
  }), names(prs))
}

test_that("a quadratic penalty's derivatives are the matrix and zero", {
  skip_if_not_installed("numDeriv")
  P <- crossprod(matrix(c(1, -1, 0, 0, 1, -1), 2, 3, byrow = TRUE))
  pen <- quadratic_penalty(P)
  b <- c(0.4, -1.2, 0.9)
  th <- list(lambda = 2.5)

  expect_equal(penalty_dhessian(pen, b, th)$lambda, unname(P),
               tolerance = 1e-12)
  expect_equal(penalty_dhessian(pen, b, th), num_dhessian(pen, b, th),
               tolerance = 1e-6)
  expect_true(all(penalty_d2hessian(pen, b, th)$lambda_lambda == 0))
  expect_true(all(penalty_dcross(pen, b, th)$lambda_lambda == 0))
  expect_true(beta_quadratic(pen, th))
  # and the Hessian really is lambda times the derivative
  expect_equal(penalty_hessian(pen, b, th),
               2.5 * penalty_dhessian(pen, b, th)$lambda, tolerance = 1e-12)
})

test_that("an additive penalty's derivatives are its components", {
  skip_if_not_installed("numDeriv")
  P1 <- crossprod(matrix(c(1, -1, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1), 3, 4,
                         byrow = TRUE))
  P2 <- diag(4)
  pen <- additive_penalty(list(P1, P2))
  b <- c(0.3, -0.7, 1.1, 0.2)
  th <- list(lambda1 = 1.7, lambda2 = 0.6)

  d <- penalty_dhessian(pen, b, th)
  expect_named(d, c("lambda1", "lambda2"))
  expect_equal(d, num_dhessian(pen, b, th), tolerance = 1e-6)
  expect_equal(penalty_hessian(pen, b, th),
               1.7 * d$lambda1 + 0.6 * d$lambda2, tolerance = 1e-12)
  expect_true(all(vapply(penalty_d2hessian(pen, b, th),
                         function(M) all(M == 0), logical(1))))
  expect_equal(penalty_dcross(pen, b, th), num_dcross(pen, b, th),
               tolerance = 1e-6)
})

test_that("a structured penalty reads the matrix parameter's derivatives", {
  skip_if_not_installed("numDeriv")
  st <- parameters7::log_cholesky(3)
  pen <- structured_penalty(st)
  b <- c(0.5, -0.4, 0.8)
  th <- as.list(stats::setNames(c(0.1, -0.2, 0.3, 0.15, -0.05, 0.2),
                                pen@params))

  expect_equal(penalty_dhessian(pen, b, th), num_dhessian(pen, b, th),
               tolerance = 1e-6)
  expect_equal(penalty_dcross(pen, b, th), num_dcross(pen, b, th),
               tolerance = 1e-6)
  # the second derivative is not zero here, the precision being a nonlinear
  # map of its free vector
  d2 <- penalty_d2hessian(pen, b, th)
  expect_true(any(vapply(d2, function(M) any(abs(M) > 1e-8), logical(1))))
  expect_true(beta_quadratic(pen, th))
})

test_that("a separable penalty carries the parent's response derivatives", {
  skip_if_not_installed("numDeriv")
  pen <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(),
                          mu = 0), n_coef = 4)
  b <- c(0.6, -1.3, 0.2, 0.9)
  th <- list(sigma = 1.4)

  d <- penalty_dhessian(pen, b, th)
  # a ridge's Hessian is I/sigma^2, so its derivative is -2 I/sigma^3, which
  # is NOT linear in sigma and is exactly what the old test of linearity used
  # to refuse
  expect_equal(d$sigma, diag(-2 / 1.4^3, 4), tolerance = 1e-10)
  expect_equal(d, num_dhessian(pen, b, th), tolerance = 1e-6)
  # exact now that the parent supplies distrib_hess_y_hess: the second
  # derivative of I/sigma^2 is 6I/sigma^4, and this used to be a difference
  expect_equal(penalty_d2hessian(pen, b, th)$sigma_sigma,
               diag(6 / 1.4^4, 4), tolerance = 1e-13)
  # and the mixed block's second derivative, 6 beta / sigma^4
  expect_equal(penalty_dcross(pen, b, th)$sigma_sigma, 6 * b / 1.4^4,
               tolerance = 1e-13)
  expect_true(beta_quadratic(pen, th))
})

test_that("a heavy-tailed prior answers too, and is not beta-quadratic", {
  skip_if_not_installed("numDeriv")
  pen <- heavy_penalty(n_coef = 3L)
  b <- c(0.7, -0.5, 1.4)
  th <- list(sigma = 1.2, nu = 5)
  expect_equal(penalty_dhessian(pen, b, th), num_dhessian(pen, b, th),
               tolerance = 1e-5)
  # a t prior's log-density is not quadratic in the coefficients
  expect_false(beta_quadratic(pen, th))
})

test_that("a kinked penalty rejects, naming what it cannot do", {
  pen <- lasso_penalty(n_coef = 3L)
  b <- c(0.2, -0.4, 0.6)
  th <- list(lambda = 1.5)
  expect_error(penalty_dhessian(pen, b, th), "has a kink")
  expect_error(penalty_d2hessian(pen, b, th), "has a kink")
  expect_error(penalty_dcross(pen, b, th), "has a kink")

  expect_error(penalty_dhessian(scad_penalty(n_coef = 3L), b,
                                list(lambda = 1, a = 3.7)),
               "does not supply")
})

test_that("the keys are the ones penalty_hess_theta uses", {
  # a consumer looks a pair up by name, so the two enumerations must agree or
  # it will read the wrong entry without anything failing
  pen <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
  b <- c(0.1, 0.2, 0.3)
  th <- list(lambda1 = 1, lambda2 = 2)
  expect_identical(names(penalty_d2hessian(pen, b, th)),
                   names(penalty_hess_theta(pen, b, th)))
  expect_identical(names(penalty_dcross(pen, b, th)),
                   names(penalty_hess_theta(pen, b, th)))
})
