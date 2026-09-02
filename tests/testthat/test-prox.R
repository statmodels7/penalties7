# The proximal operator. Every closed form is checked against the DEFINITION
# it implements -- the minimizer of (1/2t)||b - v||^2 + rho(b) -- found by a
# route that shares nothing with the formula: a direct numerical minimization
# of the objective the operator is defined by.

prox_by_search <- function(pen, v, step, theta) {
  obj <- function(b) sum((b - v)^2) / (2 * step) + penalty_value(pen, b, theta)
  # separable, so one scalar search per coordinate over a wide bracket
  vapply(seq_along(v), function(j) {
    f1 <- function(bj) {
      b <- v
      b[j] <- bj
      obj(b)
    }
    lo <- min(v[j], 0) - 10 * step - 10
    hi <- max(v[j], 0) + 10 * step + 10
    grid <- seq(lo, hi, length.out = 4001)
    b0 <- grid[which.min(vapply(grid, f1, numeric(1)))]
    stats::optimize(f1, c(b0 - (hi - lo) / 4000, b0 + (hi - lo) / 4000),
                    tol = 1e-12)$minimum
  }, numeric(1))
}

v <- c(2.3, -1.4, 0.35, 0, -0.05, 5)

test_that("the lasso prox is the soft threshold, and it is the minimizer", {
  pen <- lasso_penalty(n_coef = length(v))
  for (step in c(0.2, 1, 3)) {
    th <- list(lambda = 0.8)
    got <- penalty_prox(pen, v, step, th)
    expect_equal(got, sign(v) * pmax(abs(v) - step * 0.8, 0))
    expect_equal(got, prox_by_search(pen, v, step, th), tolerance = 1e-6)
  }
})

test_that("the ridge prox is the shrinkage, and it is the minimizer", {
  pen <- ridge_penalty(n_coef = length(v))
  # the hyperparameter is the PRECISION, so the shrinkage rises with it
  th <- list(lambda = 1 / 1.3^2)
  for (step in c(0.2, 1, 3)) {
    got <- penalty_prox(pen, v, step, th)
    expect_equal(got, v / (1 + step / 1.3^2))
    expect_equal(got, prox_by_search(pen, v, step, th), tolerance = 1e-6)
  }
})

test_that("the SCAD and MCP prox agree with the minimizer over every region", {
  vv <- c(0.1, 0.9, 1.6, 3.2, 8, -0.9, -3.2, -8)
  ps <- scad_penalty(n_coef = length(vv))
  th_s <- list(lambda = 1, a = 3.7)
  for (step in c(0.3, 1, 2)) {
    got <- penalty_prox(ps, vv, step, th_s)
    expect_equal(got, prox_by_search(ps, vv, step, th_s), tolerance = 1e-5)
  }
  # the three regions are all exercised by that grid
  expect_true(any(abs(vv) <= th_s$lambda) && any(abs(vv) > th_s$a * th_s$lambda))

  pm <- mcp_penalty(n_coef = length(vv))
  th_m <- list(lambda = 1, gamma = 3)
  for (step in c(0.3, 1, 2)) {
    got <- penalty_prox(pm, vv, step, th_m)
    expect_equal(got, prox_by_search(pm, vv, step, th_m), tolerance = 1e-5)
  }
})

test_that("a step past the concave region is rejected rather than guessed", {
  expect_error(penalty_prox(scad_penalty(n_coef = 2), c(1, 2), 2.7,
                            list(lambda = 1, a = 3.7)), "step < a - 1")
  expect_error(penalty_prox(mcp_penalty(n_coef = 2), c(1, 2), 3,
                            list(lambda = 1, gamma = 3)), "step < gamma")
})

test_that("the quadratic prox is the linear solve, with or without a map", {
  q <- 4
  P <- crossprod(diff(diag(q), differences = 1))
  pen <- quadratic_penalty(P)
  th <- list(lambda = 2)
  vq <- c(1, -2, 0.5, 3)
  got <- penalty_prox(pen, vq, 0.7, th)
  expect_equal(got, as.numeric(solve(diag(q) + 0.7 * 2 * P, vq)))
  # against the definition, minimized jointly this time
  o <- stats::optim(vq, function(b)
    sum((b - vq)^2) / (2 * 0.7) + penalty_value(pen, b, th),
    method = "BFGS", control = list(reltol = 1e-14))
  expect_equal(got, o$par, tolerance = 1e-5)

  # a map is admitted here, the objective staying quadratic
  D <- diff(diag(q))
  penD <- quadratic_penalty(diag(q - 1), map = D)
  expect_true(has_prox(penD))
  gd <- penalty_prox(penD, vq, 0.5, th)
  expect_equal(gd, as.numeric(solve(diag(q) + 0.5 * 2 * crossprod(D), vq)))
})

test_that("a separable penalty under a map is rejected, and says why", {
  D <- diff(diag(4))
  pen <- lasso_penalty(map = D)
  expect_false(has_prox(pen))
  expect_error(penalty_prox(pen, c(1, 2, 3, 4), 1, list(lambda = 1)),
               "does not split")
})

test_that("a separable penalty with no closed form is solved from its score", {
  # the heavy-tailed prior: a Student t at zero, log-concave in y for the
  # purpose of the root, and with no elementary proximal operator
  pen <- heavy_penalty(n_coef = 3)
  th <- list(sigma = 1, nu = 6)
  vt <- c(0.4, -2.5, 6)
  got <- penalty_prox(pen, vt, 0.8, th)
  # the stationary condition it claims to solve, evaluated independently
  ly <- distributions7::distrib_grad_y(pen@parent, got, th)
  expect_equal((got - vt) / 0.8, as.numeric(ly), tolerance = 1e-8)
})

test_that("a penalty without an operator refuses, and has_prox agrees", {
  st <- parameters7::log_cholesky(3)
  ps <- structured_penalty(st)
  expect_true(has_prox(ps))
  th <- stats::setNames(as.list(rep(0, length(st@free_names))), st@free_names)
  S <- penalty_hessian(ps, rep(0, 3), th)
  expect_equal(penalty_prox(ps, c(1, -1, 2), 0.5, th),
               as.numeric(solve(diag(3) + 0.5 * S, c(1, -1, 2))))
})

test_that("the generic validates its arguments", {
  pen <- lasso_penalty(n_coef = 3)
  expect_error(penalty_prox(pen, c(1, 2), 1, list(lambda = 1)), "length 3")
  expect_error(penalty_prox(pen, c(1, 2, 3), 0, list(lambda = 1)), "positive")
  expect_error(penalty_prox(pen, c(1, 2, 3), -1, list(lambda = 1)), "positive")
  expect_error(has_prox(diag(2)), "penalty object")
})


test_that("the piecewise table is the proximal operator it describes", {
  # two independent routes to the same map: the operator written out per
  # family, and the table a compiled loop reads. The table is what makes a
  # coordinate descent compilable without naming a family in the kernel, so
  # it has to agree with the operator everywhere, including across every
  # breakpoint and at the breakpoints themselves.
  q <- 4L
  cases <- list(
    ridge = list(pen = distrib_penalty(
      distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
      n_coef = q), th = list(sigma = 0.8)),
    lasso = list(pen = lasso_penalty(n_coef = q), th = list(lambda = 1.5)),
    enet  = list(pen = elasticnet_penalty(n_coef = q),
                 th = list(lambda = 1.5, alpha = 0.6)),
    scad  = list(pen = scad_penalty(n_coef = q), th = list(lambda = 1.2,
                                                           a = 3.7)),
    mcp   = list(pen = mcp_penalty(n_coef = q), th = list(lambda = 1.2,
                                                          gamma = 3))
  )
  for (nm in names(cases)) {
    pen <- cases[[nm]]$pen
    th <- cases[[nm]]$th
    for (t in c(0.05, 0.3, 1.1)) {
      sp <- penalty_prox_spec(pen, th, rep(t, q))
      expect_false(is.null(sp), label = nm)
      # a grid that crosses every breakpoint, plus the breakpoints exactly
      brk <- as.numeric(sp$cut[1L, ])
      brk <- brk[is.finite(brk)]
      u <- sort(unique(c(seq(-8, 8, by = 0.05), brk, -brk,
                         brk + 1e-12, brk - 1e-12)))
      got <- vapply(u, function(x)
        prox_apply(sp, rep(x, q))[[1L]], numeric(1))
      want <- vapply(u, function(x)
        penalty_prox(pen, rep(x, q), t, th)[[1L]], numeric(1))
      expect_equal(got, want, tolerance = 1e-12,
                   label = sprintf("%s at step %g", nm, t))
    }
  }
})

test_that("the table carries one step per coordinate", {
  # in a coordinate descent the step is 1/sum(w x^2), which differs by column
  pen <- lasso_penalty(n_coef = 3L)
  st <- c(0.2, 1, 4)
  sp <- penalty_prox_spec(pen, list(lambda = 2), st)
  expect_equal(sp$cut[, 1L], st * 2)
  for (j in seq_along(st)) {
    expect_equal(prox_apply(sp, rep(3, 3))[[j]],
                 penalty_prox(pen, rep(3, 3), st[[j]], list(lambda = 2))[[j]],
                 tolerance = 1e-12)
  }
})

test_that("the table survives a diagonal map", {
  # standardization is a diagonal map, and a penalized block reaches the
  # compiled coordinate descent through this table rather than through the
  # operator, so a standardized penalty that lost its table would silently
  # fall back on the general proximal route. Same two independent routes as
  # above, now with the map in place on both sides.
  q <- 4L
  d <- c(0.5, 0.8, 1.25, 2)
  D <- Matrix::Diagonal(x = d)
  cases <- list(
    ridge = list(pen = distrib_penalty(
      distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
      map = D), th = list(sigma = 0.8)),
    lasso = list(pen = lasso_penalty(map = D), th = list(lambda = 1.5)),
    enet  = list(pen = elasticnet_penalty(map = D),
                 th = list(lambda = 1.5, alpha = 0.6)),
    scad  = list(pen = scad_penalty(map = D), th = list(lambda = 1.2,
                                                        a = 3.7)),
    mcp   = list(pen = mcp_penalty(map = D), th = list(lambda = 1.2,
                                                       gamma = 3))
  )
  for (nm in names(cases)) {
    pen <- cases[[nm]]$pen
    th <- cases[[nm]]$th
    expect_equal(as.integer(pen@n_coef), q, label = nm)
    for (t in c(0.05, 0.2)) {
      sp <- penalty_prox_spec(pen, th, rep(t, q))
      expect_false(is.null(sp), label = nm)
      for (j in seq_len(q)) {
        brk <- as.numeric(sp$cut[j, ])
        brk <- brk[is.finite(brk)]
        u <- sort(unique(c(seq(-8, 8, by = 0.05), brk, -brk,
                           brk + 1e-12, brk - 1e-12)))
        got <- vapply(u, function(x)
          prox_apply(sp, rep(x, q))[[j]], numeric(1))
        want <- vapply(u, function(x)
          penalty_prox(pen, rep(x, q), t, th)[[j]], numeric(1))
        expect_equal(got, want, tolerance = 1e-12,
                     label = sprintf("%s at step %g, coordinate %d", nm, t, j))
      }
    }
  }
})

test_that("a diagonal map tightens the convex region of SCAD and MCP", {
  # the condition is tested on the SCALED step, t d^2 < a - 1, so a step that
  # is admissible under the identity map is not under a map that stretches
  q <- 2L
  D <- Matrix::Diagonal(x = c(1, 3))
  expect_false(is.null(penalty_prox_spec(scad_penalty(n_coef = q),
                                         list(lambda = 1, a = 3.7),
                                         rep(1, q))))
  expect_null(penalty_prox_spec(scad_penalty(map = D),
                                list(lambda = 1, a = 3.7), rep(1, q)))
  expect_false(is.null(penalty_prox_spec(mcp_penalty(n_coef = q),
                                         list(lambda = 1, gamma = 3),
                                         rep(1, q))))
  expect_null(penalty_prox_spec(mcp_penalty(map = D),
                                list(lambda = 1, gamma = 3), rep(1, q)))
})

test_that("a penalty with no such description says so", {
  # a quadratic under a general matrix is not separable, and a step past the
  # convex region of SCAD or MCP makes the operator set-valued
  expect_null(penalty_prox_spec(quadratic_penalty(diag(3)), list(lambda = 1),
                                rep(0.5, 3)))
  expect_null(penalty_prox_spec(scad_penalty(n_coef = 3L),
                                list(lambda = 1, a = 3.7), rep(3, 3)))
  expect_null(penalty_prox_spec(mcp_penalty(n_coef = 3L),
                                list(lambda = 1, gamma = 3), rep(3.5, 3)))
  # a separable penalty whose operator is a root, not a formula
  heavy <- heavy_penalty(n_coef = 3L)
  expect_null(penalty_prox_spec(heavy, list(sigma = 1, nu = 4), rep(0.5, 3)))
  # and an off-centre parent, where the map is not odd
  off <- distrib_penalty(distributions7::fixed(
    distributions7::laplace2_distrib(), mu = 0.5), n_coef = 3L, kinks = 0.5)
  expect_null(penalty_prox_spec(off, list(lambda = 1), rep(0.5, 3)))
  # a map that mixes coordinates is the generalized-lasso problem
  gen <- lasso_penalty(map = matrix(c(1, 0, 1, 1), 2, 2))
  expect_null(penalty_prox_spec(gen, list(lambda = 1), rep(0.5, 2)))
})
