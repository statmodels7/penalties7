# The step of the piecewise table, the global minimizer and the tolerance of
# the root route, and the SCAD and MCP operators near their convexity bounds.

test_that("penalty_prox_spec() rejects a step of the wrong length or value", {
  pens <- list(
    lasso = list(lasso_penalty(n_coef = 3), list(lambda = 1)),
    scad = list(scad_penalty(n_coef = 3), list(lambda = 1, a = 3.7)),
    mcp = list(mcp_penalty(n_coef = 3), list(lambda = 1, gamma = 3))
  )
  bad <- list(c(0.5, 2), numeric(0), -1, 0, c(0.5, NA), c(0.5, Inf, 1),
              rep(0.5, 4), "0.5", TRUE)
  for (p in pens) {
    for (st in bad) {
      expect_error(penalty_prox_spec(p[[1]], p[[2]], step = st),
                   "'step' must be finite and positive.*3 values")
    }
    # a single value and one value per coefficient are both accepted, and
    # agree with the operator at each coefficient's own step
    one <- penalty_prox_spec(p[[1]], p[[2]], step = 0.5)
    three <- penalty_prox_spec(p[[1]], p[[2]], step = c(0.5, 0.25, 1))
    expect_equal(dim(one$cut), dim(three$cut))
    u <- c(1.7, -2.4, 0.9)
    for (j in 1:3) {
      st <- c(0.5, 0.25, 1)[j]
      expect_equal(prox_apply(three, u)[j],
                   penalty_prox(p[[1]], u, st, p[[2]])[j], tolerance = 1e-14)
    }
  }
})

# the global minimizer of (b - v)^2/(2t) - log f(b) for a Student t of scale
# s and nu degrees of freedom, found without the package: every local minimum
# on a fine grid (a change of sign of the closed-form derivative from - to +),
# polished by uniroot, and the one of smallest objective kept
t_obj <- function(x, v, t, s, nu) (x - v)^2 / (2 * t) - log(dt(x / s, nu) / s)
t_dobj <- function(x, v, t, s, nu) (x - v) / t + (nu + 1) * x / (nu * s^2 + x^2)
t_prox_ref <- function(v, t, s, nu) {
  g <- sort(unique(c(seq(min(0, v) - 1, max(0, v) + 1, length.out = 100001),
                     seq(min(0, v), max(0, v), length.out = 50001), 0)))
  d <- t_dobj(g, v, t, s, nu)
  k <- which(d[-length(d)] < 0 & d[-1] > 0)
  cand <- c(g[d == 0], vapply(k, function(i)
    stats::uniroot(t_dobj, c(g[i], g[i + 1]), v = v, t = t, s = s, nu = nu,
                   tol = 1e-300, maxiter = 5000)$root, numeric(1)))
  cand[which.min(t_obj(cand, v, t, s, nu))]
}

test_that("the root route returns the global minimizer for a Student t prior", {
  h <- heavy_penalty(n_coef = 1)
  # the two measured cases: the stationary point the bracket found first was
  # a local minimum of larger objective
  expect_equal(penalty_prox(h, 12, 8, list(sigma = 1, nu = 4)), 2,
               tolerance = 1e-6)
  r2 <- t_prox_ref(3, 1, 0.3, 1)
  expect_equal(r2, 0.16698253579062719, tolerance = 1e-9)
  expect_equal(penalty_prox(h, 3, 1, list(sigma = 0.3, nu = 1)), r2,
               tolerance = 1e-6)

  # cases from a scan of 384, the first six of which returned the other
  # local minimum before
  cs <- rbind(c(-3, 1, 0.3, 1), c(12, 8, 0.3, 2), c(-12, 8, 0.3, 2),
              c(40, 40, 0.3, 4), c(-40, 40, 0.3, 4), c(-12, 8, 1, 4),
              c(0.5, 40, 0.3, 4), c(-40, 0.1, 3, 30), c(12, 1, 1, 1))
  for (i in seq_len(nrow(cs))) {
    th <- list(sigma = cs[i, 3], nu = cs[i, 4])
    expect_equal(penalty_prox(h, cs[i, 1], cs[i, 2], th),
                 t_prox_ref(cs[i, 1], cs[i, 2], cs[i, 3], cs[i, 4]),
                 tolerance = 1e-8)
  }
  # the same points as one vector, each coordinate solved on its own
  pv <- heavy_penalty(n_coef = 3)
  th <- list(sigma = 0.3, nu = 2)
  v <- c(12, -0.5, 3)
  expect_equal(penalty_prox(pv, v, 8, th),
               vapply(v, t_prox_ref, numeric(1), t = 8, s = 0.3, nu = 2),
               tolerance = 1e-8)
})

test_that("the root route is accurate relative to a tiny answer", {
  h <- heavy_penalty(n_coef = 1)
  th <- list(sigma = 1, nu = 4)
  # in the linear regime l_y(b) = -(nu + 1) b / (nu sigma^2), so the answer
  # is v / (1 + t (nu + 1) / (nu sigma^2)) to a relative O(v^2)
  for (v in c(1e-11, -1e-20, 1e-200)) {
    expect_lt(abs(penalty_prox(h, v, 1, th) / (v / (1 + 5 / 4)) - 1), 1e-10)
  }
  h3 <- heavy_penalty(n_coef = 3)
  v <- c(1e-11, 1e-20, -1e-200)
  expect_lt(max(abs(penalty_prox(h3, v, 1, th) / (v / (1 + 5 / 4)) - 1)), 1e-10)
})

test_that("the SCAD operator is exact as the step approaches a - 1", {
  sc <- scad_penalty(n_coef = 1)
  th <- list(lambda = 1, a = 3.7)
  # step = 2.7 (1 - 1e-12) and u = 3.69999999999865 as doubles; the exact
  # value of the operator at these doubles, from mpmath at 60 digits, rounds
  # to the double 2.35
  st <- (3.7 - 1) * (1 - 1e-12)
  u <- 3.69999999999865
  expect_identical(sprintf("%a", st), "0x1.59999999981dap+1")
  expect_identical(sprintf("%a", u), "0x1.d999999998dbap+1")
  expect_equal(penalty_prox(sc, u, st, th), 2.35, tolerance = 1e-13)
  expect_equal(penalty_prox(sc, -u, st, th), -2.35, tolerance = 1e-13)
  # the curvature-scaled operator, at step st/2 and curvature 2, is the
  # unscaled one at step st and rate lambda / 2
  sc2 <- scad_penalty(n_coef = 1, curvature = 2)
  expect_equal(penalty_prox(sc2, u, st / 2, list(lambda = 2, a = 3.7)), 2.35,
               tolerance = 1e-13)

  # away from the bound the middle piece agrees with the textbook form
  uu <- c(1.6, 2.2, 3.1, -2.9, 3.7)
  for (step in c(0.3, 1, 2)) {
    textbook <- sign(uu) * (abs(uu) - step * 3.7 / 2.7) / (1 - step / 2.7)
    mid <- abs(uu) > 1 + step & abs(uu) <= 3.7
    got <- penalty_prox(scad_penalty(n_coef = 5), uu, step, th)
    expect_equal(got[mid], textbook[mid], tolerance = 4 * .Machine$double.eps)
  }
})

test_that("the MCP operator is exact as the step approaches gamma", {
  mc <- mcp_penalty(n_coef = 1)
  th <- list(lambda = 1, gamma = 3)
  # step = 3 (1 - 1e-12) and u = 3 - 1e-13 as doubles; the exact value at
  # these doubles, from mpmath at 60 digits, is 2.900074019245003701
  st <- 3 * (1 - 1e-12)
  u <- 3 - 1e-13
  expect_identical(sprintf("%a", st), "0x1.7ffffffffe59dp+1")
  expect_identical(sprintf("%a", u), "0x1.7ffffffffff1fp+1")
  expect_equal(penalty_prox(mc, u, st, th), 2.900074019245003701,
               tolerance = 1e-13)
  # away from the bound it agrees with the textbook form
  uu <- c(0.5, 1.2, 2.2, -2.9, 3)
  for (step in c(0.3, 1, 2)) {
    textbook <- sign(uu) * pmax(abs(uu) - step, 0) / (1 - step / 3)
    expect_equal(penalty_prox(mcp_penalty(n_coef = 5), uu, step, th), textbook,
                 tolerance = 4 * .Machine$double.eps)
  }
})
