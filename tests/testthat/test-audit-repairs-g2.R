# The matrix of a quadratic penalty must be positive semidefinite, and the
# rank and the null basis are read from the same eigenvalues.

test_that("quadratic_penalty() rejects a P that is not positive semidefinite", {
  # diag(2, -1) was built with rank 1, and its value at (0, 100) was
  # 1/2 * (-10000) + constant, unbounded below
  expect_error(quadratic_penalty(diag(c(2, -1))),
               "'P' is not positive semidefinite.", fixed = TRUE)
  expect_error(quadratic_penalty(-diag(2)),
               "'P' is zero or negative semidefinite.", fixed = TRUE)
  expect_error(quadratic_penalty(matrix(0, 3, 3)),
               "'P' is zero or negative semidefinite.", fixed = TRUE)
  expect_error(quadratic_penalty(diag(c(2, -1)), map = diag(2)),
               "'P' is not positive semidefinite.", fixed = TRUE)
  expect_error(quadratic_penalty(diag(c(2, -1)), blocks = 3),
               "'P' is not positive semidefinite.", fixed = TRUE)
  expect_error(quadratic_penalty(-diag(2), blocks = 3),
               "'P' is zero or negative semidefinite.", fixed = TRUE)

  # a rank-deficient P is accepted, and so is one whose null eigenvalues
  # come out of the arithmetic slightly negative
  pen <- quadratic_penalty(crossprod(diff(diag(5), differences = 2)))
  expect_identical(penalty_rank(pen), 3L)
  set.seed(3)
  Q <- qr.Q(qr(matrix(rnorm(16), 4)))
  P <- Q %*% diag(c(1, 0.5, 0, -1e-14)) %*% t(Q)
  P <- (P + t(P)) / 2
  expect_lt(min(eigen(P, symmetric = TRUE, only.values = TRUE)$values), 0)
  pen <- quadratic_penalty(P)
  expect_identical(penalty_rank(pen), 2L)
  expect_identical(ncol(penalty_null_basis(pen)), 2L)
})

test_that("quadratic_penalty() rejects a map that annihilates P", {
  # D'PD = 0 exactly; the penalty was built with rank 0 and value 0
  P <- diag(c(0, 0, 1))
  D <- rbind(c(1, 0), c(0, 1), c(0, 0))
  expect_error(quadratic_penalty(P, map = D),
               "'map' annihilates 'P'", fixed = TRUE)
  expect_error(quadratic_penalty(P, map = Matrix::Matrix(D, sparse = TRUE)),
               "'map' annihilates 'P'", fixed = TRUE)
  # a map that leaves part of P is accepted
  D2 <- rbind(c(1, 0), c(0, 1), c(1, 1))
  expect_identical(penalty_rank(quadratic_penalty(P, map = D2)), 1L)
})

test_that("the rank and the null basis of a quadratic penalty agree", {
  # an eigenvalue at the threshold of the rank rule: with two decompositions
  # 68 of these 3000 penalties had q - rank != ncol(null basis)
  set.seed(11)
  mism <- 0L
  for (i in seq_len(3000)) {
    Q <- qr.Q(qr(matrix(rnorm(16), 4)))
    P <- Q %*% diag(c(1, 0.3, 1e-10 * (1 + 1e-6 * rnorm(1)), 0)) %*% t(Q)
    P <- (P + t(P)) / 2
    pen <- quadratic_penalty(P)
    if (4L - penalty_rank(pen) != ncol(penalty_null_basis(pen))) {
      mism <- mism + 1L
    }
  }
  expect_identical(mism, 0L)

  # without a map the rank and the constant are the ones of the values-only
  # decomposition of P, and the null basis is orthonormal and annihilated
  set.seed(4)
  A <- matrix(rnorm(36), 6)
  mats <- list(crossprod(diff(diag(6))),
               crossprod(diff(diag(6), differences = 2)),
               crossprod(A),
               crossprod(A[1:3, ]))
  for (P in mats) {
    P <- (P + t(P)) / 2
    pen <- quadratic_penalty(P)
    ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
    keep <- ev > 1e-10 * max(ev, 0)
    expect_identical(pen@logpdet_DPD, sum(log(ev[keep])))
    expect_identical(penalty_rank(pen), sum(keep))
    expect_identical(penalty_logpdet(pen, list(lambda = 2))$value,
                     sum(keep) * log(2) + sum(log(ev[keep])))
    N <- penalty_null_basis(pen)
    expect_identical(ncol(N), 6L - sum(keep))
    expect_equal(crossprod(N), diag(ncol(N)), tolerance = 1e-12)
    expect_lt(max(abs(P %*% N), 0), 1e-12 * max(ev))
  }
})
