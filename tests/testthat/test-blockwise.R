# A distrib_penalty whose parent is multivariate is read one BLOCK at a time:
# the coordinates of a block depend on each other and the blocks do not.

.mv_centered <- function(p, ...) {
  do.call(distributions7::fixed,
          c(list(distributions7::mvgaussian_distrib(p, ...)),
            stats::setNames(as.list(rep(0, p)), paste0("mu", seq_len(p)))))
}


test_that("the block width comes from the parent and the count must divide", {
  pen <- distrib_penalty(.mv_centered(2), n_coef = 6)
  expect_identical(pen@block, 2L)
  expect_identical(pen@n_coef, 6L)
  expect_error(distrib_penalty(.mv_centered(2), n_coef = 7), "whole blocks")

  # a univariate parent is the same construction at a block of one, so the
  # route that existed before is unchanged rather than special-cased
  u <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 5)
  expect_identical(u@block, 1L)
})


test_that("the value is the negative log density of the blocks", {
  # the reference is written out from the multivariate normal density and
  # shares no arithmetic with the penalty
  set.seed(7)
  for (p in 2:3) {
    m <- 4
    pen <- distrib_penalty(.mv_centered(p), n_coef = m * p)
    eta <- stats::rnorm(length(pen@params), sd = 0.3)
    th <- stats::setNames(as.list(eta), pen@params)
    b <- stats::rnorm(m * p)
    sig <- distributions7::mv_sigma(
      distributions7::mvgaussian_distrib(p),
      c(stats::setNames(as.list(rep(0, p)), paste0("mu", seq_len(p))), th))
    B <- matrix(b, ncol = p, byrow = TRUE)
    ref <- sum(apply(B, 1, function(r) {
      p / 2 * log(2 * pi) +
        determinant(sig, logarithm = TRUE)$modulus[1] / 2 +
        sum(r * solve(sig, r)) / 2
    }))
    expect_equal(penalty_value(pen, b, th), ref, tolerance = 1e-12,
                 info = paste("p =", p))
  }
})


test_that("the Hessian is block diagonal, exactly", {
  # the blocks are independent, so a non-zero entry across two of them would
  # be a coupling the model does not have
  pen <- distrib_penalty(.mv_centered(2), n_coef = 8)
  th <- list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1, sigma_L2.1 = 0.5)
  b <- c(0.3, -0.4, 1.1, 0.2, -0.9, 0.6, 0.1, -0.2)
  h <- penalty_hessian(pen, b, th)
  for (i in 1:4) {
    ix <- (i - 1L) * 2L + 1:2
    h[ix, ix] <- 0
  }
  expect_equal(max(abs(h)), 0)
})


test_that("every derivative of a blockwise penalty agrees with numDeriv", {
  skip_if_not_installed("numDeriv")
  set.seed(7)
  for (p in 2:3) {
    m <- 4
    pen <- distrib_penalty(.mv_centered(p), n_coef = m * p)
    nm <- pen@params
    eta <- stats::rnorm(length(nm), sd = 0.3)
    th <- stats::setNames(as.list(eta), nm)
    b <- stats::rnorm(m * p)

    fb <- function(z) penalty_value(pen, z, th)
    expect_equal(penalty_gradient(pen, b, th), numDeriv::grad(fb, b),
                 tolerance = 1e-6, info = paste("p =", p))
    expect_equal(unname(penalty_hessian(pen, b, th)), numDeriv::hessian(fb, b),
                 tolerance = 1e-4)

    ft <- function(z) penalty_value(pen, b, stats::setNames(as.list(z), nm))
    expect_equal(unname(unlist(penalty_grad_theta(pen, b, th))),
                 numDeriv::grad(ft, eta), tolerance = 1e-6)
    hn <- numDeriv::hessian(ft, eta)
    hh <- penalty_hess_theta(pen, b, th)
    for (i in seq_along(nm)) for (j in seq_along(nm)) {
      key <- if (i == j) paste0(nm[i], "_", nm[i]) else
        paste(nm[sort(c(i, j))], collapse = "_")
      if (!is.null(hh[[key]])) {
        expect_equal(hh[[key]], hn[i, j], tolerance = 1e-4,
                     info = paste("p =", p, key))
      }
    }
    fg <- function(z) penalty_gradient(pen, b, stats::setNames(as.list(z), nm))
    expect_equal(unname(do.call(cbind, penalty_cross(pen, b, th))),
                 numDeriv::jacobian(fg, eta), tolerance = 1e-6)
  }
})


test_that("a block whose covariance is diagonal IS a product of univariates", {
  # the twin: independent coordinates within a block make the blockwise
  # reading and the coordinatewise one the same prior, so the two routes agree
  # with no tolerance to choose. The chart of diagonal_matrix is the log of the
  # DIAGONAL ENTRY, i.e. of the variance, so the univariate scale is its root.
  mv <- do.call(distributions7::fixed,
    list(distributions7::mvgaussian_distrib(
      2, sigma = parameters7::diagonal_matrix(2)), mu1 = 0, mu2 = 0))
  pb <- distrib_penalty(mv, n_coef = 6)
  th <- stats::setNames(as.list(c(log(1.5), log(0.7))), pb@params)
  b <- c(0.3, -0.8, 1.2, 0.4, -0.5, 0.9)

  u <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 3)
  ref <- penalty_value(u, b[c(1, 3, 5)], list(sigma = sqrt(1.5))) +
         penalty_value(u, b[c(2, 4, 6)], list(sigma = sqrt(0.7)))
  expect_equal(penalty_value(pb, b, th), ref, tolerance = 1e-12)
})


test_that("a blockwise penalty has no proximal operator and says so", {
  # the operator acts one coordinate at a time and a correlated block does not
  # separate; the predicate is what a fitting layer routes on
  pen <- distrib_penalty(.mv_centered(2), n_coef = 6)
  th <- list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.3)
  expect_false(has_prox(pen))
  expect_null(penalty_prox_spec(pen, th, 1))
  expect_error(penalty_prox(pen, rep(0, 6), 1, th), "do not separate")
  # and a kink is a point of a scalar argument, so there is none
  expect_length(penalty_kinks(pen, th), 0L)
})


test_that("a blockwise penalty passes check_penalty", {
  skip_if_not_installed("numDeriv")
  pen <- distrib_penalty(.mv_centered(3), n_coef = 9)
  th <- stats::setNames(as.list(c(0.1, -0.2, 0.15, 0.3, -0.1, 0.2)),
                        pen@params)
  res <- check_penalty(pen, theta = th, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})


test_that("a blockwise penalty says what its hyperparameters are ABOUT", {
  # the hyperparameters are the free values of a matrix parameter -- the
  # logarithms of a Cholesky diagonal and the entries below it -- and nobody
  # reads those; the quantities are the standard deviations and correlations
  pen <- distrib_penalty(.mv_centered(2), n_coef = 6)
  th <- list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1, sigma_L2.1 = 0.5)
  rd <- penalty_readable(pen, th)
  expect_identical(names(rd$value), c("sd_v1", "sd_v2", "cor_v1_v2"))
  expect_identical(unname(rd$transform), c("log", "log", "atanh"))
  expect_identical(dim(rd$jacobian), c(3L, 3L))

  # the quantities ARE those of the prior's own family
  expect_equal(rd$value,
               distributions7::mv_derived(pen@parent, th)$value)

  # and the Jacobian is right: one central difference on the quantities
  f <- function(v) {
    distributions7::mv_derived(
      pen@parent, stats::setNames(as.list(v), pen@params))$value
  }
  v0 <- unlist(th)
  num <- vapply(seq_along(v0), function(k) {
    up <- dn <- v0
    up[k] <- v0[k] + 1e-5
    dn[k] <- v0[k] - 1e-5
    (f(up) - f(dn)) / 2e-5
  }, numeric(3))
  expect_equal(unname(rd$jacobian), unname(num), tolerance = 1e-6)
})


test_that("a penalty whose hyperparameters ARE the quantities answers NULL", {
  # that is the honest answer for every other branch: a smoothing parameter, a
  # rate, a shape are each read on their own scale already
  expect_null(penalty_readable(quadratic_penalty(diag(2)), list(lambda = 1)))
  expect_null(penalty_readable(lasso_penalty(n_coef = 2), list(lambda = 1)))
  expect_null(penalty_readable(scad_penalty(n_coef = 2),
                               list(lambda = 1, a = 3.7)))
  u <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    n_coef = 3)
  expect_null(penalty_readable(u, list(sigma = 1)))
})


test_that("the marginal pieces of a blockwise penalty agree with numDeriv", {
  skip_if_not_installed("numDeriv")
  # what reml() and ml() read to estimate a correlated random effect's
  # covariance: the theta derivative of the penalty Hessian, its second, and
  # the theta Hessian of the mixed block. Each is one difference of the
  # ANALYTIC quantity below it.
  set.seed(13)
  for (p in 2:3) {
    m <- 3
    pen <- distrib_penalty(.mv_centered(p), n_coef = m * p)
    nm <- pen@params
    eta <- stats::rnorm(length(nm), sd = 0.3)
    th <- stats::setNames(as.list(eta), nm)
    b <- stats::rnorm(m * p)
    at <- function(z) stats::setNames(as.list(z), nm)

    dh <- penalty_dhessian(pen, b, th)
    expect_named(dh, nm)
    for (k in seq_along(nm)) {
      ref <- matrix(numDeriv::jacobian(function(z) {
        vv <- eta; vv[k] <- z
        as.vector(penalty_hessian(pen, b, at(vv)))
      }, eta[k]), m * p, m * p)
      expect_equal(dh[[k]], ref, tolerance = 1e-5,
                   info = paste("dhessian p =", p, nm[k]))
    }

    d2h <- penalty_d2hessian(pen, b, th)
    dc <- penalty_dcross(pen, b, th)
    prs <- names(d2h)
    for (nmk in prs) {
      ij <- strsplit(nmk, "_")[[1L]]
      a <- match(ij[1L], nm)
      bb <- match(ij[length(ij)], nm)
      if (is.na(a) || is.na(bb)) next
      ref2 <- matrix(numDeriv::jacobian(function(z) {
        vv <- eta; vv[bb] <- z
        as.vector(penalty_dhessian(pen, b, at(vv))[[a]])
      }, eta[bb]), m * p, m * p)
      expect_equal(d2h[[nmk]], ref2, tolerance = 1e-4,
                   info = paste("d2hessian p =", p, nmk))
      ref3 <- numDeriv::jacobian(function(z) {
        vv <- eta; vv[bb] <- z
        penalty_cross(pen, b, at(vv))[[a]]
      }, eta[bb])[, 1L]
      expect_equal(dc[[nmk]], ref3, tolerance = 1e-4,
                   info = paste("dcross p =", p, nmk))
    }
  }
})


test_that("the marginal pieces work for a Student t parent too", {
  skip_if_not_installed("numDeriv")
  # the same three quantities, on a family whose response Hessian DEPENDS on
  # the observation: the blockwise assembly reads one matrix per block there
  # rather than one matrix, and dp_blockdiag carries both shapes
  set.seed(31)
  p <- 2
  m <- 3
  mvt <- do.call(distributions7::fixed,
                 list(distributions7::mvstudent_t_distrib(p),
                      mu1 = 0, mu2 = 0))
  pen <- distrib_penalty(mvt, n_coef = m * p)
  nm <- pen@params
  expect_identical(nm, c("sigma_log_L1", "sigma_log_L2", "sigma_L2.1", "nu"))
  eta <- c(stats::rnorm(3, sd = 0.3), 7)
  th <- stats::setNames(as.list(eta), nm)
  b <- stats::rnorm(m * p)
  at <- function(z) stats::setNames(as.list(z), nm)

  dh <- penalty_dhessian(pen, b, th)
  for (k in seq_along(nm)) {
    ref <- matrix(numDeriv::jacobian(function(z) {
      vv <- eta; vv[k] <- z
      as.vector(penalty_hessian(pen, b, at(vv)))
    }, eta[k]), m * p, m * p)
    expect_equal(dh[[k]], ref, tolerance = 1e-5, info = nm[k])
  }
  d2h <- penalty_d2hessian(pen, b, th)
  dc <- penalty_dcross(pen, b, th)
  for (nmk in names(d2h)) {
    ij <- strsplit(nmk, "_")[[1L]]
    a <- match(ij[1L], nm)
    bb <- match(ij[length(ij)], nm)
    if (is.na(a) || is.na(bb)) next
    ref2 <- matrix(numDeriv::jacobian(function(z) {
      vv <- eta; vv[bb] <- z
      as.vector(penalty_dhessian(pen, b, at(vv))[[a]])
    }, eta[bb]), m * p, m * p)
    expect_equal(d2h[[nmk]], ref2, tolerance = 1e-4, info = nmk)
    ref3 <- numDeriv::jacobian(function(z) {
      vv <- eta; vv[bb] <- z
      penalty_cross(pen, b, at(vv))[[a]]
    }, eta[bb])[, 1L]
    expect_equal(dc[[nmk]], ref3, tolerance = 1e-4, info = nmk)
  }
})
