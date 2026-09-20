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
test_that("the additive branch is quadratic and answers as one", {
  # A sum of quadratic forms is a quadratic form. is_quadratic() inherited
  # FALSE from the base class while penalty_matrix(), penalty_rank() and
  # penalty_logpdet() all answered and penalty_null_basis() rejected, pointing
  # the reader at a predicate that would have told them nothing.
  P1 <- crossprod(diff(diag(5)))
  P2 <- crossprod(diff(diag(5), differences = 2))
  pen <- additive_penalty(list(P1, P2))

  expect_true(is_quadratic(pen))

  # the null space of the sum is the intersection of the components', which for
  # these two is the constants, and it does not move with the hyperparameters
  nb <- penalty_null_basis(pen)
  expect_identical(dim(nb), c(5L, 1L))
  for (th in list(list(lambda1 = 1, lambda2 = 1),
                  list(lambda1 = 1e-6, lambda2 = 1e6))) {
    S <- penalty_matrix(pen, th)
    expect_lt(max(abs(S %*% nb)) / max(abs(S)), 1e-12)
  }
  # the count is an integer, sum() over a logical, as it was before this
  expect_identical(penalty_rank(pen), 4L)

  # the log pseudo-determinant answers in the shape the other two quadratic
  # branches use: named lists keyed by hyperparameter and by pair. It returned
  # an unnamed numeric vector and a matrix, so the row that now reaches it,
  # which reads lp$grad[[p]], would have been an error rather than a number.
  th <- list(lambda1 = 0.7, lambda2 = 2.5)
  lp <- penalty_logpdet(pen, th)
  expect_identical(names(lp$grad), pen@params)
  expect_identical(names(lp$hess), names(penalty_hess_theta(
    pen, rep(0.1, pen@n_coef), th)))
  expect_type(lp$grad, "list")
  expect_type(lp$hess, "list")

  for (p in pen@params) {
    ref <- numDeriv::grad(function(v) {
      t2 <- th; t2[[p]] <- v
      penalty_logpdet(pen, t2)$value
    }, th[[p]])
    expect_equal(lp$grad[[p]], ref, tolerance = 1e-7)
  }

  # check_penalty() runs the three quadratic rows on the branch now, and they
  # pass: 8 rows before, 11 after
  res <- check_penalty(pen, theta = th, verbose = FALSE)
  expect_identical(nrow(res), 11L)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
  expect_true(all(c("quadratic three-point identity",
                    "logpdet gradient vs numDeriv",
                    "null basis annihilates the matrix") %in% res$check))

  # and the proximal route is unchanged: has_prox() asks whether penalty_prox
  # is registered on the base class before it asks is_quadratic(), and this
  # branch registers none
  expect_false(has_prox(pen))
})
test_that("a decomposition that cannot resolve its smallest kept eigenvalue is rejected", {
  P1 <- crossprod(diff(diag(5)))                    # first differences,  rank 4
  P2 <- crossprod(diff(diag(5), differences = 2))   # second differences, rank 3
  b <- c(0.4, -1.1, 0.7, 0.2, -0.3)

  # The components are ordered so that the one carrying lambda1 is rank
  # deficient ON THE RANGE: the four range directions are covered by the
  # second differences only three at a time, so the fourth is carried by
  # lambda2 alone and does not grow with lambda1. Written the other way round
  # the spread produces no small eigenvalue at all and nothing degenerates.
  pen <- additive_penalty(list(P2, P1))
  r <- penalty_rank(pen)
  expect_identical(r, 4L)

  # An ordinary setting is untouched, and is the eigen route written out.
  th <- list(lambda1 = 2, lambda2 = 0.5)
  S <- penalty_matrix(pen, th)
  e <- eigen(S, symmetric = TRUE)
  k <- order(e$values, decreasing = TRUE)[seq_len(r)]
  expect_equal(penalty_logpdet(pen, th)$value, sum(log(e$values[k])))
  expect_true(all(is.finite(unlist(penalty_grad_theta(pen, b, th)))))

  # A wide spread is not on its own a degenerate one. At 1e10 the smallest
  # range eigenvalue is still four thousand times the resolution and the
  # answer is finite, which is what keeps the rejection from reading the
  # spread in place of the conditioning.
  wide <- list(lambda1 = 1e10, lambda2 = 1)
  ew <- sort(eigen(penalty_matrix(pen, wide), symmetric = TRUE,
                   only.values = TRUE)$values, decreasing = TRUE)
  expect_gt(ew[r], 1e3 * .Machine$double.eps * ew[1] * pen@n_coef)
  expect_true(is.finite(penalty_logpdet(pen, wide)$value))

  # Past the precision of a double it is. The smallest range eigenvalue is
  # carried by lambda2 alone, so it stays at ew[r] while the resolution grows
  # with lambda1: the premise is a gap of sixteen orders and not a sign, so
  # the verdict does not turn on platform arithmetic.
  far <- list(lambda1 = 1e30, lambda2 = 1)
  res_far <- .Machine$double.eps *
    max(eigen(penalty_matrix(pen, far), symmetric = TRUE,
              only.values = TRUE)$values) * pen@n_coef
  expect_gt(res_far / ew[r], 1e10)

  # Everything the decomposition feeds is NaN there, and no warning is
  # raised: a search visiting such a point would otherwise raise one per
  # evaluation.
  expect_silent(v <- penalty_value(pen, b, far))
  expect_true(is.nan(v))
  lp <- penalty_logpdet(pen, far)
  expect_true(is.nan(lp$value))
  expect_true(all(is.nan(unlist(lp$grad))))
  expect_true(all(is.nan(unlist(lp$hess))))
  expect_true(all(is.nan(unlist(penalty_grad_theta(pen, b, far)))))
  expect_true(all(is.nan(unlist(penalty_hess_theta(pen, b, far)))))

  # The pseudo-inverse is rejected whole rather than left as it came out,
  # which is what an injection has to fail: without the rejection it is a
  # finite matrix and every derivative above is a finite number.
  a <- penalties7:::additive_sum(pen, far)
  expect_true(all(is.nan(a$Sp)))
  expect_identical(dim(a$Sp), c(pen@n_coef, pen@n_coef))

  # The sum itself is assembled without a decomposition, so it is returned as
  # it stands and the coefficient derivatives stay finite.
  expect_true(all(is.finite(a$S)))
  expect_true(all(is.finite(penalty_gradient(pen, b, far))))
  expect_true(all(is.finite(penalty_hessian(pen, b, far))))

  # The negative control: a large penalty is not a degenerate one. Raising
  # both parameters together leaves the condition number where it was, so
  # nothing is rejected and the value moves by the scaling alone.
  big <- list(lambda1 = 1e30, lambda2 = 1e30)
  one <- list(lambda1 = 1, lambda2 = 1)
  expect_true(is.finite(penalty_logpdet(pen, big)$value))
  expect_equal(penalty_logpdet(pen, big)$value - penalty_logpdet(pen, one)$value,
               r * log(1e30))
})
