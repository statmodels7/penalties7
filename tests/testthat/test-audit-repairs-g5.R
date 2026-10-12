# Repairs of the additive penalty: a component annihilated by the map, and
# settings at which the rank identity cannot be evaluated.

test_that("a component that the map carries to zero is rejected by number", {
  P <- list(diag(3), diag(c(0, 0, 1)))
  D0 <- rbind(c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0), c(0, 0, 0, 0, 0))
  expect_error(additive_penalty(P, map = D0),
               "component 2 is zero after 'map' is applied")

  # the same components through a map that keeps the third row
  D <- rbind(c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
  pen <- additive_penalty(P, map = D)
  expect_identical(pen@n_coef, 5L)
  expect_identical(pen@p_rank, 3L)
  expect_equal(pen@mats[[2]], t(D) %*% P[[2]] %*% D)
})

test_that("a setting whose rank identity is NA returns NaN without an error", {
  P1 <- crossprod(diff(diag(5)))
  P2 <- crossprod(diff(diag(5), differences = 2))
  pen <- additive_penalty(list(P1, P2))
  b <- c(0.4, -1.1, 0.7, 0.2, -0.3)

  for (l1 in c(1e-308, 1e-310)) {
    th <- list(lambda1 = l1, lambda2 = 1)
    expect_no_error(v <- penalty_value(pen, b, th))
    expect_true(is.nan(v))
    expect_no_error(g <- penalty_grad_theta(pen, b, th))
    expect_true(all(is.nan(unlist(g))))
    expect_no_error(h <- penalty_hess_theta(pen, b, th))
    expect_true(all(is.nan(unlist(h))))
    expect_no_error(lp <- penalty_logpdet(pen, th))
    expect_true(all(is.nan(unlist(lp))))
  }

  # At 1e-306 the point is still resolved. To first order in lambda1 the log
  # pseudo-determinant is that of P2 plus log(lambda1 u'P1u), with u the unit
  # linear trend: the direction P2 leaves free and P1 does not.
  l1 <- 1e-306
  th <- list(lambda1 = l1, lambda2 = 1)
  u <- (1:5 - 3) / sqrt(10)
  ev2 <- eigen(P2, symmetric = TRUE, only.values = TRUE)$values[1:3]
  ref <- sum(log(ev2)) + log(l1) + log(sum(u * (P1 %*% u)))
  lp <- penalty_logpdet(pen, th)
  expect_equal(lp$value, ref, tolerance = 1e-12)
  expect_equal(l1 * lp$grad$lambda1, 1, tolerance = 1e-8)
  expect_true(is.finite(penalty_value(pen, b, th)))
  expect_true(all(is.finite(unlist(penalty_grad_theta(pen, b, th)))))
})

test_that("a weighted sum or a reduction that overflows returns NaN", {
  P1 <- crossprod(diff(diag(5)))
  P2 <- crossprod(diff(diag(5), differences = 2))
  pen <- additive_penalty(list(P1, P2))
  b <- c(0.4, -1.1, 0.7, 0.2, -0.3)

  # 1e308 * P2 overflows the sum itself; at the second setting the sum is
  # finite and the reduction on the null space of the dominant group is not
  for (th in list(list(lambda1 = 1, lambda2 = 1e308),
                  list(lambda1 = 5e-324, lambda2 = 1e-310))) {
    expect_no_error(lp <- penalty_logpdet(pen, th))
    expect_true(all(is.nan(unlist(lp))))
    expect_no_error(v <- penalty_value(pen, b, th))
    expect_true(is.nan(v))
    expect_no_error(g <- penalty_grad_theta(pen, b, th))
    expect_true(all(is.nan(unlist(g))))
  }
})
