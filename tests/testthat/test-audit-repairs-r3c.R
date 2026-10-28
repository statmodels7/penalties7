# The rows of check_penalty() that compare the stored rank with the
# eigenvalues of the matrix and the width of the null basis with the rank.

r3c_dd <- function(q, k) crossprod(diff(diag(q), differences = k))

r3c_with <- function(pen, ...) {
  do.call(QuadraticPenalty, utils::modifyList(S7::props(pen), list(...)))
}

r3c_rows <- c("rank vs eigenvalues of the matrix",
              "null basis width vs n_coef minus rank")

test_that("check_penalty() fails a stored rank one below or one above the true one", {
  skip_if_not_installed("numDeriv")
  # second differences over five coefficients: rank 3, two null directions
  good <- quadratic_penalty(r3c_dd(5, 2))
  expect_identical(penalty_rank(good), 3L)
  for (r in c(2, 4)) {
    bad <- r3c_with(good, rank = r)
    expect_identical(penalty_rank(bad), as.integer(r))
    res <- check_penalty(bad, verbose = FALSE)
    for (nm in r3c_rows) {
      row <- res[res$check == nm, ]
      expect_identical(nrow(row), 1L, label = paste(r, nm))
      expect_identical(row$status, "FAILED", label = paste(r, nm))
      expect_identical(row$max_error, 1, label = paste(r, nm))
    }
    # every other row reads the stored rank and passes
    expect_true(all(res$status[!res$check %in% r3c_rows] == "OK"),
                label = paste("other rows at rank", r))
  }
})

test_that("check_penalty() fails a null basis narrower than n_coef minus the rank", {
  skip_if_not_installed("numDeriv")
  good <- quadratic_penalty(r3c_dd(5, 2))
  bad <- r3c_with(good, null_basis = good@null_basis[, 1L, drop = FALSE])
  res <- check_penalty(bad, verbose = FALSE)
  st <- stats::setNames(res$status, res$check)
  expect_identical(unname(st[r3c_rows]), c("OK", "FAILED"))
  expect_identical(res$max_error[res$check == r3c_rows[2L]], 1)
  expect_identical(unname(st["null basis annihilates the matrix"]), "OK")
})

test_that("every shipped quadratic branch passes the rank rows at the defaults", {
  skip_if_not_installed("numDeriv")
  # each penalty with its rank from the closed form: a difference matrix of
  # order k over q coefficients has rank q - k
  cases <- list(
    ridge = list(ridge_penalty(n_coef = 4), 4L),
    d1 = list(quadratic_penalty(r3c_dd(5, 1)), 4L),
    d2 = list(quadratic_penalty(r3c_dd(5, 2)), 3L),
    d3_q40 = list(quadratic_penalty(r3c_dd(40, 3)), 37L),
    map = list(quadratic_penalty(diag(5), map = diff(diag(6))), 5L),
    map_ridge = list(ridge_penalty(map = diff(diag(6), differences = 2)), 4L),
    blocks = list(quadratic_penalty(r3c_dd(5, 2), blocks = 3), 9L),
    ar1 = list(structured_penalty(parameters7::ar1(4)), 4L),
    compound_symmetry = list(
      structured_penalty(parameters7::compound_symmetry(3)), 3L),
    log_cholesky = list(structured_penalty(parameters7::log_cholesky(3)), 3L),
    covariance = list(
      structured_penalty(parameters7::inverse_of(parameters7::ar1(4))), 4L),
    scaled_d2 = list(
      structured_penalty(parameters7::scaled_matrix(r3c_dd(5, 2))), 3L),
    additive2_full = list(additive_penalty(list(diag(c(1, 1, 0)),
                                                diag(c(0, 1, 1)))), 3L),
    additive2_null = list(additive_penalty(list(r3c_dd(5, 1),
                                                r3c_dd(5, 2))), 4L),
    additive3_full = list(additive_penalty(list(r3c_dd(5, 1), diag(5),
                                                r3c_dd(5, 2))), 5L),
    additive3_null = list(additive_penalty(list(r3c_dd(6, 1), r3c_dd(6, 2),
                                                r3c_dd(6, 3))), 5L)
  )
  for (nm in names(cases)) {
    pen <- cases[[nm]][[1L]]
    expect_identical(penalty_rank(pen), cases[[nm]][[2L]], label = nm)
    res <- check_penalty(pen, verbose = FALSE)
    expect_true(all(r3c_rows %in% res$check), label = nm)
    expect_identical(res$max_error[res$check %in% r3c_rows], c(0, 0),
                     label = nm)
    expect_true(all(res$status == "OK"),
                label = paste(nm, paste(res$check[res$status != "OK"],
                                        collapse = ", ")))
  }
})

test_that("check_penalty() has the row counts stated on its page", {
  skip_if_not_installed("numDeriv")
  n <- function(pen) nrow(check_penalty(pen, verbose = FALSE))
  expect_identical(n(lasso_penalty(n_coef = 4)), 5L)
  expect_identical(n(scad_penalty(n_coef = 4)), 9L)
  expect_identical(n(ridge_penalty(n_coef = 4)), 9L)
  expect_identical(n(quadratic_penalty(r3c_dd(5, 2))), 10L)
  expect_identical(
    n(structured_penalty(parameters7::log_cholesky(3))),
    39L)
  expect_identical(n(additive_penalty(list(r3c_dd(5, 1), r3c_dd(5, 2)))), 14L)
  expect_identical(n(additive_penalty(list(diag(c(1, 1, 0)),
                                           diag(c(0, 1, 1))))), 13L)
})

test_that("the rank row counts at the tolerance stored in the penalty", {
  skip_if_not_installed("numDeriv")
  # one eigenvalue at 1e-8 relative: a zero at tol = 1e-6, not at 1e-10
  set.seed(3)
  Q <- qr.Q(qr(matrix(stats::rnorm(16), 4)))
  quad <- quadratic_penalty(Q %*% diag(c(3, 1, 3e-8, 0)) %*% t(Q), tol = 1e-6)
  add <- additive_penalty(list(diag(c(1, 0, 0)), diag(c(0, 1, 1e-8))),
                          tol = 1e-6)
  expect_identical(c(quad@tol, add@tol), c(1e-6, 1e-6))
  expect_identical(c(penalty_rank(quad), penalty_rank(add)), c(2L, 2L))
  for (pen in list(quad, add)) {
    res <- check_penalty(pen, verbose = FALSE)
    expect_identical(res$status[res$check %in% r3c_rows], c("OK", "OK"))
    expect_true(all(res$status == "OK"))
  }
  # the blocked path stores the tolerance too, and the default is unchanged
  expect_identical(quadratic_penalty(r3c_dd(5, 2), blocks = 3, tol = 1e-6)@tol,
                   1e-6)
  expect_identical(quadratic_penalty(r3c_dd(5, 2))@tol, 1e-10)
})

test_that("the rank row of an additive penalty reads the normalized components", {
  skip_if_not_installed("numDeriv")
  # components whose scales differ by more than the tolerance: S(theta) has
  # one eigenvalue above 1e-10 relative, the family has rank 2
  for (s in c(1e-11, 1e-12, 1e-14)) {
    pen <- additive_penalty(list(diag(c(1, 0, 0)), diag(c(0, s, 0))))
    expect_identical(penalty_rank(pen), 2L)
    res <- check_penalty(pen, verbose = FALSE)
    row <- res[res$check == r3c_rows[1L], ]
    expect_identical(row$status, "OK", label = paste("scale", s))
    expect_identical(row$max_error, 0, label = paste("scale", s))
  }
  # a stored rank one below the count still fails both rows
  good <- additive_penalty(list(r3c_dd(5, 1), r3c_dd(5, 2)))
  bad <- do.call(AdditivePenalty,
                 utils::modifyList(S7::props(good), list(rank = 3)))
  res <- check_penalty(bad, verbose = FALSE)
  for (nm in r3c_rows) {
    row <- res[res$check == nm, ]
    expect_identical(row$status, "FAILED", label = nm)
    expect_identical(row$max_error, 1, label = nm)
  }
})

test_that("a stored rank above the count gives a table, not an error", {
  skip_if_not_installed("numDeriv")
  # the log pseudo-determinant takes the log of a zero eigenvalue, so the
  # value is NaN and numDeriv stops on every reference that reads it
  good <- additive_penalty(list(r3c_dd(5, 1), r3c_dd(5, 2)))
  bad <- do.call(AdditivePenalty,
                 utils::modifyList(S7::props(good), list(rank = 5)))
  res <- check_penalty(bad, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 14L)
  st <- stats::setNames(res$status, res$check)
  er <- stats::setNames(res$max_error, res$check)
  expect_identical(unname(st[r3c_rows]), c("FAILED", "FAILED"))
  expect_identical(unname(er[r3c_rows]), c(1, 1))
  nan_rows <- c("gradient vs numDeriv", "grad_theta[lambda1] vs numDeriv",
                "grad_theta[lambda2] vs numDeriv",
                "logpdet gradient vs numDeriv")
  expect_true(all(is.nan(er[nan_rows])))
  expect_true(all(st[is.nan(er)] == "FAILED"))
  # the rows whose reference does not read the value are still computed
  expect_identical(unname(st["hessian vs numDeriv on the gradient"]), "OK")
})

test_that("a reference that stops fails its row, an analytic error is signaled", {
  skip_if_not_installed("numDeriv")
  # a value that is not finite away from the given beta: numDeriv stops on
  # the gradient row's reference
  NaOff <- S7::new_class("NaOff", parent = QuadraticPenalty)
  b0 <- c(0.4, -0.7, 1.1)
  na <- do.call(NaOff, S7::props(quadratic_penalty(diag(3))))
  S7::method(penalty_value, NaOff) <- function(penalty, beta, theta, ...) {
    if (!identical(beta, b0)) return(NA_real_)
    S7::method(penalty_value, QuadraticPenalty)(penalty, beta, theta, ...)
  }
  res <- check_penalty(na, beta = b0, verbose = FALSE)
  row <- res[res$check == "gradient vs numDeriv", ]
  expect_true(is.nan(row$max_error))
  expect_identical(row$status, "FAILED")
  expect_identical(res$status[res$check == "hessian vs numDeriv on the gradient"],
                   "OK")
  # an analytic method that stops at beta and theta still signals its error
  ErrGrad <- S7::new_class("ErrGrad", parent = QuadraticPenalty)
  eg <- do.call(ErrGrad, S7::props(quadratic_penalty(diag(3))))
  S7::method(penalty_gradient, ErrGrad) <- function(penalty, beta, theta, ...) {
    stop("broken gradient")
  }
  expect_error(check_penalty(eg, verbose = FALSE), "broken gradient")
  ErrValue <- S7::new_class("ErrValue", parent = QuadraticPenalty)
  ev <- do.call(ErrValue, S7::props(quadratic_penalty(diag(3))))
  S7::method(penalty_value, ErrValue) <- function(penalty, beta, theta, ...) {
    stop("broken value")
  }
  expect_error(check_penalty(ev, verbose = FALSE), "broken value")
  expect_error(check_penalty(quadratic_penalty(diag(3)),
                             theta = list(lambda = -1)), "open interval")
})
