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

# beta_quadratic compares the parent's response Hessian BETWEEN readings and
# not between the entries of one of them, and it probes at four POINTS rather
# than at four values. Each half has its own negative control below, written
# out as the form it replaces, so neither can be satisfied by weakening the
# predicate.

# penalties7 0.21.0: every entry of one reading against the first entry
bq_entrywise <- function(pen, theta) {
  t <- c(-1.73, -0.29, 0.61, 2.04)
  h <- tryCatch(distributions7::distrib_hess_y(pen@parent, t,
                                               align_ptheta(pen, theta)),
                error = function(e) NULL)
  !is.null(h) && length(h) > 1L &&
    isTRUE(all.equal(as.numeric(h), rep(as.numeric(h)[1L], length(h)),
                     tolerance = 1e-12))
}

# the same comparison between readings, but probed at four VALUES laid out
# with pen@block columns rather than at four points
bq_four_values <- function(pen, theta) {
  h <- tryCatch(distributions7::distrib_hess_y(
    pen@parent, suppressWarnings(dp_arg(pen, c(-1.73, -0.29, 0.61, 2.04))),
    align_ptheta(pen, theta)), error = function(e) NULL)
  if (is.null(h) || !length(h)) return(FALSE)
  k <- as.integer(pen@block)^2L
  if (length(h) %% k) return(FALSE)
  n <- length(h) %/% k
  if (n <= 1L) return(k > 1L || length(h) > 1L)
  v <- as.numeric(h)
  isTRUE(all.equal(v, rep(v[seq_len(k)], n), tolerance = 1e-12))
}

# the mean is flattened into mu1..mup, so fixed() is given one name each
mv_zero_mean <- function(d, p) {
  do.call(distributions7::fixed,
          c(list(d), as.list(stats::setNames(rep(0, p),
                                             paste0("mu", seq_len(p))))))
}
mv_theta <- function(d) {
  as.list(stats::setNames(rep(0, length(d@params)), d@params))
}

test_that("a multivariate gaussian prior is beta-quadratic at every width", {
  for (p in 2:4) {
    par <- mv_zero_mean(distributions7::mvgaussian1_distrib(
      p, parameters7::log_cholesky(p)), p)
    pen <- distrib_penalty(par, n_coef = 4L * p)
    th <- mv_theta(par)
    expect_identical(pen@block, p)
    # hess_y is the constant -Sigma^{-1}, so the prior is exactly quadratic
    expect_true(beta_quadratic(pen, th))
    # and the entrywise form compares an off-diagonal with a diagonal, so it
    # answers FALSE precisely where the prior is quadratic
    expect_false(bq_entrywise(pen, th))
  }
})

test_that("a multivariate t prior is not, at a width four values cannot see", {
  pens <- lapply(c(2L, 4L), function(p) {
    par <- mv_zero_mean(distributions7::mvstudent_t1_distrib(
      p, parameters7::log_cholesky(p)), p)
    th <- mv_theta(par)
    th[["nu"]] <- 6
    list(pen = distrib_penalty(par, n_coef = 4L * p), th = th)
  })
  for (cs in pens) {
    # hess_y carries the (nu + p)/(nu + q) reweighting, so it moves with the
    # response and the prior is quadratic at no width
    expect_false(beta_quadratic(cs$pen, cs$th))
  }
  # the negative control here is the WIDTH. Four values laid out with p
  # columns are two points at p = 2 and ONE from p = 4 up, and a single
  # reading takes the branch that answers without comparing anything: the
  # four-value probe sees the p = 2 case and reports the p = 4 one quadratic.
  expect_false(bq_four_values(pens[[1L]]$pen, pens[[1L]]$th))
  expect_true(bq_four_values(pens[[2L]]$pen, pens[[2L]]$th))
})

test_that("a univariate parent's answer did not move", {
  g <- distrib_penalty(distributions7::fixed(
    distributions7::gaussian1_distrib(), mu = 0), n_coef = 4L)
  cases <- list(list(g, list(sigma = 1.3), TRUE),
                list(heavy_penalty(n_coef = 4L), list(sigma = 1, nu = 4),
                     FALSE),
                list(lasso_penalty(n_coef = 4L), list(lambda = 1), TRUE))
  for (cs in cases) {
    expect_identical(cs[[1L]]@block, 1L)
    expect_identical(beta_quadratic(cs[[1L]], cs[[2L]]), cs[[3L]])
    # a block of one coordinate reads the same under both forms, which is
    # what says the change is confined to the multivariate branch
    expect_identical(bq_entrywise(cs[[1L]], cs[[2L]]), cs[[3L]])
  }
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

# penalty_dhessian_beta() ----------------------------------------------------
#
# The reference differentiates penalty_hessian() along the direction with ONE
# central difference, so it shares no arithmetic with the closed form, which
# reads the parent's third response derivative.

num_dhessian_beta <- function(pen, beta, theta, v, h = 1e-5) {
  (penalty_hessian(pen, beta + h * v, theta) -
     penalty_hessian(pen, beta - h * v, theta)) / (2 * h)
}

test_that("a heavy-tailed prior's Hessian moves with the coefficients", {
  pen <- heavy_penalty(n_coef = 4L)
  b <- c(0.8, -1.9, 0.3, 2.6)
  v <- c(0.4, 0.7, -1.1, 0.2)
  th <- list(sigma = 0.9, nu = 3.5)
  got <- penalty_dhessian_beta(pen, b, th, v)
  ref <- num_dhessian_beta(pen, b, th, v)
  expect_true(is.matrix(got))
  expect_identical(dim(got), c(4L, 4L))
  expect_gt(max(abs(ref)), 1e-2)
  expect_equal(got, unname(as.matrix(ref)), tolerance = 1e-7)
  # linear in the direction, which a contraction must be
  expect_equal(penalty_dhessian_beta(pen, b, th, 2 * v), 2 * got,
               tolerance = 1e-12)
})

test_that("a map is carried on both sides and on the direction", {
  D <- rbind(c(1, -1, 0), c(0, 1, -1), c(0.5, 0, 1), c(0, 2, 0))
  pen <- distrib_penalty(
    distributions7::fixed(distributions7::student_t1_distrib(), mu = 0),
    map = D)
  b <- c(0.3, -0.8, 1.4)
  v <- c(-0.5, 0.9, 0.25)
  th <- list(sigma = 1.3, nu = 5)
  expect_equal(penalty_dhessian_beta(pen, b, th, v),
               unname(as.matrix(num_dhessian_beta(pen, b, th, v))),
               tolerance = 1e-7)
})

test_that("a penalty quadratic in the coefficients answers zero", {
  b <- c(0.2, -0.4, 0.6)
  v <- c(1, 2, 3)
  z <- matrix(0, 3, 3)
  expect_identical(penalty_dhessian_beta(quadratic_penalty(diag(3)), b,
                                         list(lambda = 2), v), z)
  expect_identical(penalty_dhessian_beta(
    additive_penalty(list(diag(3), diag(c(1, 1, 0)))), b,
    list(lambda1 = 1, lambda2 = 2), v), z)
  s <- structured_penalty(parameters7::log_cholesky(3))
  expect_identical(penalty_dhessian_beta(
    s, b, as.list(stats::setNames(rep(0.1, 6), s@params)), v), z)
  g <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 3L)
  expect_identical(penalty_dhessian_beta(g, b, list(sigma = 1.4), v), z)
})

test_that("a kinked or a non-quadratic multivariate parent rejects", {
  b <- c(0.2, -0.4, 0.6)
  v <- c(1, 0, -1)
  expect_error(penalty_dhessian_beta(lasso_penalty(n_coef = 3L), b,
                                     list(lambda = 1), v), "has a kink")
  expect_error(penalty_dhessian_beta(scad_penalty(n_coef = 3L), b,
                                     list(lambda = 1, a = 3.7), v),
               "does not supply")
  expect_error(penalty_dhessian_beta(heavy_penalty(n_coef = 3L), b,
                                     list(sigma = 1, nu = 4), c(1, 2)),
               "has length 2")
})
