# SCAD and MCP scaled by the curvature of the loss in each coordinate:
# sum_j c_j p(t_j; lambda / c_j). The slope at zero stays lambda and the knee
# moves to a lambda / c_j.

one_coord <- function(mk, t, lam, extra) {
  pen <- mk(n_coef = 1L)
  penalty_value(pen, t, c(list(lambda = lam), extra))
}

test_that("the scaled value is the sum of unit penalties at lambda / c", {
  t <- c(-2.3, -0.4, 0.05, 0.9, 3.9, 12)
  cv <- c(0.5, 2, 3.7, 10, 40, 1)
  lam <- 1.3
  for (case in list(list(mk = scad_penalty, extra = list(a = 3.7)),
                    list(mk = mcp_penalty, extra = list(gamma = 2.6)))) {
    pen <- case$mk(n_coef = 6L, curv = cv)
    th <- c(list(lambda = lam), case$extra)
    ref <- sum(vapply(seq_along(t), function(j)
      cv[j] * one_coord(case$mk, t[j], lam / cv[j], case$extra), numeric(1)))
    expect_equal(penalty_value(pen, t, th), ref, tolerance = 1e-12)
    # the slope at zero is lambda in every coordinate
    g0 <- penalty_gradient(pen, rep(1e-9, 6), th)
    expect_equal(g0, rep(lam, 6), tolerance = 1e-6)
  }
})

test_that("unit curvature is the unscaled penalty to the last bit", {
  t <- c(-2.3, -0.4, 0.05, 0.9, 3.9)
  for (case in list(list(mk = scad_penalty, extra = list(a = 3.7)),
                    list(mk = mcp_penalty, extra = list(gamma = 2.6)))) {
    th <- c(list(lambda = 1.3), case$extra)
    a <- case$mk(n_coef = 5L)
    b <- case$mk(n_coef = 5L, curv = 1)
    expect_identical(penalty_value(a, t, th), penalty_value(b, t, th))
    expect_identical(penalty_gradient(a, t, th), penalty_gradient(b, t, th))
    expect_identical(penalty_grad_theta(a, t, th), penalty_grad_theta(b, t, th))
  }
})

test_that("check_penalty passes with a curvature", {
  cv <- c(0.5, 2, 3.7, 10)
  expect_true(all(check_penalty(scad_penalty(n_coef = 4L, curv = cv),
                                theta = list(lambda = 0.8, a = 3.7))$ok))
  expect_true(all(check_penalty(mcp_penalty(n_coef = 4L, curv = cv),
                                theta = list(lambda = 0.8, gamma = 2.6))$ok))
})

test_that("the operator and its table follow the curvature", {
  v <- seq(-6, 6, length.out = 41)
  cv <- rep(c(0.5, 2, 3.7), length.out = length(v))
  for (case in list(list(mk = scad_penalty, extra = list(a = 3.7)),
                    list(mk = mcp_penalty, extra = list(gamma = 2.6)))) {
    th <- c(list(lambda = 1.3), case$extra)
    pen <- case$mk(n_coef = length(v), curv = cv)
    step <- rep(0.4, length(v))
    got <- penalty_prox(pen, v, 0.4, th)
    # coordinatewise, the unscaled operator at step s c and rate lambda / c
    ref <- vapply(seq_along(v), function(j)
      penalty_prox(case$mk(n_coef = 1L), v[j], step[j] * cv[j],
                   c(list(lambda = 1.3 / cv[j]), case$extra)), numeric(1))
    expect_equal(got, ref, tolerance = 1e-12)
    sp <- penalty_prox_spec(pen, th, step)
    expect_false(is.null(sp))
    # reading the table reproduces the operator
    tab <- vapply(seq_along(v), function(j) {
      k <- which(abs(v[j]) <= sp$cut[j, ])[1L]
      sign(v[j]) * (sp$slope[j, k] * abs(v[j]) + sp$icept[j, k])
    }, numeric(1))
    expect_equal(tab, got, tolerance = 1e-12)
  }
})

test_that("a step of 1 / c is admissible at any shape the bounds allow", {
  cv <- c(0.01, 1, 100)
  pen <- scad_penalty(n_coef = 3L, curv = cv)
  expect_false(is.null(penalty_prox_spec(pen, list(lambda = 1, a = 2.05),
                                         1 / cv)))
  pen <- mcp_penalty(n_coef = 3L, curv = cv)
  expect_false(is.null(penalty_prox_spec(pen, list(lambda = 1, gamma = 1.05),
                                         1 / cv)))
  expect_error(scad_penalty(n_coef = 3L, curv = c(1, -1, 2)), "positive")
  expect_error(scad_penalty(n_coef = 3L, curv = c(1, 2)), "3 values")
})
