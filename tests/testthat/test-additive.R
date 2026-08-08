# The sum of quadratic penalties, one smoothing parameter per component.

two_way <- function(k = 4) {
  D <- crossprod(diff(diag(k), differences = 2))
  list(kronecker(diag(k), D), kronecker(D, diag(k)))
}

test_that("the value and the derivatives are what the definition says", {
  skip_if_not_installed("numDeriv")
  mats <- two_way(4)
  pen <- additive_penalty(mats)
  expect_identical(pen@params, c("lambda1", "lambda2"))

  set.seed(2)
  beta <- rnorm(16)
  th <- list(lambda1 = 0.7, lambda2 = 12)

  # the value, written out here from the definition
  S <- th$lambda1 * mats[[1]] + th$lambda2 * mats[[2]]
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  r <- penalty_rank(pen)
  lpd <- sum(log(sort(ev, decreasing = TRUE)[seq_len(r)]))
  expect_equal(penalty_value(pen, beta, th),
               0.5 * sum(beta * (S %*% beta)) - 0.5 * lpd + r / 2 * log(2 * pi))

  # every derivative against a numerical one of the value
  expect_equal(penalty_gradient(pen, beta, th),
               numDeriv::grad(function(b) penalty_value(pen, b, th), beta),
               tolerance = 1e-6)
  expect_equal(penalty_hessian(pen, beta, th),
               numDeriv::hessian(function(b) penalty_value(pen, b, th), beta),
               tolerance = 1e-5)

  gt <- unlist(penalty_grad_theta(pen, beta, th))
  num_gt <- numDeriv::grad(function(l)
    penalty_value(pen, beta, list(lambda1 = l[1], lambda2 = l[2])),
    c(th$lambda1, th$lambda2))
  expect_equal(unname(gt), num_gt, tolerance = 1e-6)

  ht <- penalty_hess_theta(pen, beta, th)
  num_ht <- numDeriv::hessian(function(l)
    penalty_value(pen, beta, list(lambda1 = l[1], lambda2 = l[2])),
    c(th$lambda1, th$lambda2))
  expect_equal(ht[["lambda1_lambda1"]], num_ht[1, 1], tolerance = 1e-4)
  expect_equal(ht[["lambda2_lambda2"]], num_ht[2, 2], tolerance = 1e-4)
  expect_equal(ht[["lambda1_lambda2"]], num_ht[1, 2], tolerance = 1e-4)

  cr <- penalty_cross(pen, beta, th)
  num_cr <- numDeriv::jacobian(function(b)
    unlist(penalty_grad_theta(pen, b, th)), beta)
  expect_equal(cr[["lambda1"]], num_cr[1, ], tolerance = 1e-6)
  expect_equal(cr[["lambda2"]], num_cr[2, ], tolerance = 1e-6)
})

test_that("the rank does not move as the parameters spread apart", {
  # The measured trap: counting eigenvalues of the assembled sum reads the
  # rank correctly only while the components are comparable, and loses the
  # small contributions once they are not. The stored rank is a property of
  # the components and cannot move.
  mats <- two_way(4)
  pen <- additive_penalty(mats)
  r <- penalty_rank(pen)

  counted <- vapply(c(1, 1e4, 1e10, 1e14), function(ratio) {
    S <- mats[[1]] + ratio * mats[[2]]
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    sum(ev > 1e-10 * max(ev))
  }, numeric(1))
  expect_true(any(counted < r))          # the count really does fall
  expect_identical(penalty_rank(pen), r) # and the object's rank does not

  # the null basis annihilates the sum at every parameter value, which is
  # the statement the rank stands for
  stacked <- Reduce(`+`, lapply(mats, function(P) P / max(abs(P))))
  e <- eigen(stacked, symmetric = TRUE)
  N <- e$vectors[, e$values <= 1e-10 * max(e$values), drop = FALSE]
  expect_gt(ncol(N), 0)
  for (lam in list(c(1, 1), c(1e-6, 1e6))) {
    S <- lam[1] * mats[[1]] + lam[2] * mats[[2]]
    expect_lt(max(abs(S %*% N)) / max(abs(S)), 1e-12)
  }
})

test_that("anisotropy is what the sum buys over one scaled matrix", {
  mats <- two_way(4)
  pen <- additive_penalty(mats)
  set.seed(5)
  beta <- rnorm(16)
  # penalizing one direction hard and the other not is a value a single
  # scaled matrix cannot produce, whatever its scale
  v_aniso <- penalty_value(pen, beta, list(lambda1 = 1e-3, lambda2 = 1e3))
  v_iso <- penalty_value(pen, beta, list(lambda1 = 1, lambda2 = 1))
  expect_false(isTRUE(all.equal(v_aniso, v_iso)))
  # and the quadratic form differs in the direction that is smoothed
  expect_gt(sum(beta * (mats[[2]] %*% beta)) * 1e3,
            sum(beta * (mats[[1]] %*% beta)) * 1e-3)
})

test_that("the whole contract passes check_penalty", {
  pen <- additive_penalty(two_way(4))
  res <- check_penalty(pen, beta = rnorm(16),
                       theta = list(lambda1 = 0.5, lambda2 = 3),
                       verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a map is carried into the components, and inputs are validated", {
  D <- diff(diag(6))
  pen <- additive_penalty(list(diag(5), crossprod(diff(diag(5)))), map = D)
  expect_identical(pen@n_coef, 6L)
  expect_identical(length(pen@params), 2L)

  expect_error(additive_penalty(list()), "non-empty")
  expect_error(additive_penalty(list(diag(3), diag(4))), "same dimension")
  expect_error(additive_penalty(list(matrix(c(1, 2, 3, 4), 2))), "symmetric")
  expect_error(additive_penalty(list(-diag(3))), "negative definite")
})

test_that("a single component reproduces the quadratic penalty", {
  P <- crossprod(diff(diag(6), differences = 2))
  a <- additive_penalty(list(P))
  q <- quadratic_penalty(P)
  set.seed(7)
  b <- rnorm(6)
  expect_equal(penalty_value(a, b, list(lambda1 = 2)),
               penalty_value(q, b, list(lambda = 2)))
  expect_equal(penalty_rank(a), penalty_rank(q))
})
