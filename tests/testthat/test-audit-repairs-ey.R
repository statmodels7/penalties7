# The tolerance of an additive penalty reaches the similarity transformation,
# the message for a matrix without a positive eigenvalue names the
# semidefinite case, and the constructors reject entries that are not finite.

test_that("additive_penalty() passes its tol to the rank count of additive_stable()", {
  # S(1e8, 1) = diag(1e8 + 1, 1e-3), rank 2 at tol = 1e-12. Reference values
  # from mpmath at 50 digits on the same double entries.
  pen <- additive_penalty(list(diag(c(1, 1e-11)), diag(c(1, 0))), tol = 1e-12)
  expect_identical(pen@tol, 1e-12)
  expect_identical(penalty_rank(pen), 2L)

  lp <- function(l1, l2) penalty_logpdet(pen, list(lambda1 = l1, lambda2 = l2))
  at <- lp(1e8, 1)
  expect_equal(at$value, 11.512925474970228, tolerance = 1e-12)
  expect_equal(at$grad$lambda1 / 1.9999999900000001e-8, 1, tolerance = 1e-12)
  expect_equal(at$grad$lambda2 / 9.9999999000000010e-9, 1, tolerance = 1e-12)

  # central differences: the value for the gradient, the gradient for the
  # Hessian; the step in lambda2 stays below lambda2 = 1
  h1 <- 1e4
  h2 <- 0.5
  cd_g1 <- (lp(1e8 + h1, 1)$value - lp(1e8 - h1, 1)$value) / (2 * h1)
  cd_g2 <- (lp(1e8, 1 + h2)$value - lp(1e8, 1 - h2)$value) / (2 * h2)
  expect_equal(cd_g1 / at$grad$lambda1, 1, tolerance = 1e-6)
  expect_equal(cd_g2 / at$grad$lambda2, 1, tolerance = 1e-5)
  cd_h11 <- (lp(1e8 + h1, 1)$grad$lambda1 - lp(1e8 - h1, 1)$grad$lambda1) /
    (2 * h1)
  cd_h22 <- (lp(1e8, 1 + h2)$grad$lambda2 - lp(1e8, 1 - h2)$grad$lambda2) /
    (2 * h2)
  cd_h12 <- (lp(1e8, 1 + h2)$grad$lambda1 - lp(1e8, 1 - h2)$grad$lambda1) /
    (2 * h2)
  expect_equal(cd_h11 / at$hess$lambda1_lambda1, 1, tolerance = 1e-6)
  expect_equal(cd_h22 / at$hess$lambda2_lambda2, 1, tolerance = 1e-6)
  expect_equal(cd_h12 / at$hess$lambda1_lambda2, 1, tolerance = 1e-6)

  # the default stores 1e-10, and an explicit 1e-10 gives the same results
  # to the last bit
  P1 <- crossprod(diff(diag(5)))
  P2 <- crossprod(diff(diag(5), differences = 2))
  a0 <- additive_penalty(list(P1, P2, diag(c(1, 0, 0, 0, 0))))
  a1 <- additive_penalty(list(P1, P2, diag(c(1, 0, 0, 0, 0))), tol = 1e-10)
  expect_identical(a0@tol, 1e-10)
  b <- c(0.4, -1.1, 0.7, 0.2, -0.3)
  for (th in list(list(lambda1 = 2, lambda2 = 0.5, lambda3 = 1),
                  list(lambda1 = 1, lambda2 = 1e9, lambda3 = 1e-6))) {
    expect_identical(penalty_logpdet(a0, th), penalty_logpdet(a1, th))
    expect_identical(penalty_value(a0, b, th), penalty_value(a1, b, th))
    expect_identical(penalty_hess_theta(a0, b, th),
                     penalty_hess_theta(a1, b, th))
  }
})

test_that("a negative semidefinite matrix is named as such by both constructors", {
  N <- diag(c(0, -1, 0))
  expect_error(quadratic_penalty(N),
               "'P' is zero or negative semidefinite.", fixed = TRUE)
  expect_error(quadratic_penalty(N, blocks = 2),
               "'P' is zero or negative semidefinite.", fixed = TRUE)
  expect_error(additive_penalty(list(diag(3), N)),
               "component 2 is zero or negative semidefinite.", fixed = TRUE)
})

test_that("entries that are not finite are rejected with the argument named", {
  A <- diag(3)
  A[2, 2] <- NA
  B <- diag(3)
  B[1, 1] <- Inf
  expect_error(additive_penalty(list(diag(3), A)),
               "component 2 of 'matrices' must have only finite entries.",
               fixed = TRUE)
  expect_error(additive_penalty(list(B, diag(3))),
               "component 1 of 'matrices' must have only finite entries.",
               fixed = TRUE)
  expect_error(quadratic_penalty(A),
               "'P' must have only finite entries.", fixed = TRUE)
  expect_error(quadratic_penalty(B),
               "'P' must have only finite entries.", fixed = TRUE)
  expect_error(quadratic_penalty(A, blocks = 2),
               "'P' must have only finite entries.", fixed = TRUE)

  M <- diag(3)
  M[1, 2] <- NA
  Mi <- diag(3)
  Mi[1, 2] <- -Inf
  for (D in list(M, Mi, Matrix::Diagonal(3, c(1, NaN, 2)),
                 Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(1, Inf, 1)))) {
    expect_error(quadratic_penalty(diag(3), map = D),
                 "'map' must have only finite entries.", fixed = TRUE)
    expect_error(additive_penalty(list(diag(3)), map = D),
                 "'map' must have only finite entries.", fixed = TRUE)
  }

  # a unit diagonal Matrix stores no entries and is accepted
  expect_identical(quadratic_penalty(diag(3), map = Matrix::Diagonal(3))@rank,
                   3L)
  expect_identical(additive_penalty(list(diag(3)),
                                    map = Matrix::Diagonal(3))@rank, 3L)
})
