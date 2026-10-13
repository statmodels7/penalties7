# The root route of the separable proximal operator for a parent whose mode
# is not at zero.

# the global minimizer of (b - v)^2/(2t) - log f(b) for a Student t of
# location mu, scale s and nu degrees of freedom, found without the package:
# a fine grid over a range holding v and mu, denser near mu and near v, every
# local minimum of the grid polished by optimize(), the smallest one kept
tm_obj <- function(b, v, t, mu, s, nu) {
  (b - v)^2 / (2 * t) - log(dt((b - mu) / s, nu) / s)
}
tm_prox_ref <- function(v, t, mu, s, nu) {
  g <- sort(unique(c(seq(min(v, mu) - 1, max(v, mu) + 1, length.out = 200001),
                     mu + s * seq(-60, 60, length.out = 24001),
                     v + seq(-1, 1, length.out = 2001))))
  f <- tm_obj(g, v, t, mu, s, nu)
  k <- unique(c(which(diff(sign(diff(f))) > 0) + 1L, which.min(f)))
  cand <- vapply(k, function(i) {
    stats::optimize(tm_obj, c(g[max(i - 1L, 1L)], g[min(i + 1L, length(g))]),
                    v = v, t = t, mu = mu, s = s, nu = nu, tol = 1e-15)$minimum
  }, numeric(1))
  cand[which.min(tm_obj(cand, v, t, mu, s, nu))]
}
t1_at <- function(mu) {
  distrib_penalty(distributions7::fixed(distributions7::student_t1_distrib(),
                                        mu = mu), n_coef = 1)
}

test_that("the root route returns the global minimizer for a parent not centered at zero", {
  p <- t1_at(-3)
  # the three measured cases: the scan from 0 to v returned the stationary
  # point near 1 (or sqrt(5)), whose objective is larger by about 8 to 12
  cs <- rbind(c(0.001, 4), c(0.01, 4), c(0.001, 2))
  for (i in seq_len(nrow(cs))) {
    th <- list(sigma = cs[i, 1], nu = 1)
    got <- penalty_prox(p, 3, cs[i, 2], th)
    ref <- tm_prox_ref(3, cs[i, 2], -3, cs[i, 1], 1)
    expect_lt(abs(got - ref), 1e-8)
    expect_lt(tm_obj(got, 3, cs[i, 2], -3, cs[i, 1], 1) -
                tm_obj(ref, 3, cs[i, 2], -3, cs[i, 1], 1), 1e-9)
  }
  expect_equal(penalty_prox(p, 3, 4, list(sigma = 0.001, nu = 1)), -2.99999925,
               tolerance = 1e-9)

  # random cases over the location, the scale, the degrees of freedom, the
  # step and the point
  set.seed(11)
  for (i in 1:25) {
    mu <- stats::runif(1, -10, 10)
    s <- exp(stats::runif(1, log(0.001), log(10)))
    nu <- exp(stats::runif(1, log(0.5), log(30)))
    t <- exp(stats::runif(1, log(0.01), log(100)))
    v <- stats::runif(1, -40, 40)
    got <- penalty_prox(t1_at(mu), v, t, list(sigma = s, nu = nu))
    ref <- tm_prox_ref(v, t, mu, s, nu)
    fr <- tm_obj(ref, v, t, mu, s, nu)
    expect_lt(tm_obj(got, v, t, mu, s, nu) - fr, 1e-9 * max(1, abs(fr)))
  }

  # the coordinates of one vector share the mode
  p3 <- distrib_penalty(
    distributions7::fixed(distributions7::student_t1_distrib(), mu = -3),
    n_coef = 3)
  v <- c(3, -10, -2.5)
  expect_equal(penalty_prox(p3, v, 4, list(sigma = 0.01, nu = 1)),
               vapply(v, tm_prox_ref, numeric(1), t = 4, mu = -3, s = 0.01,
                      nu = 1),
               tolerance = 1e-8)
})

test_that("a centered parent has its mode at zero exactly and the scan from 0 to v", {
  # the scan from 0 to v, as the route ran it for every parent, with the
  # widened bracket where the products of the scan underflow to zero
  scan0 <- function(pen, v, step, theta, n_grid = 64L) {
    theta <- align_ptheta(pen, theta)
    ly <- function(b) {
      as.numeric(distributions7::distrib_grad_y(pen@parent, b, theta)) + 0 * b
    }
    q <- length(v)
    G <- outer(seq(0, 1, length.out = n_grid), v)
    H <- (G - rep(v, each = n_grid)) / step -
      matrix(ly(as.numeric(G)), nrow = n_grid, ncol = q)
    vapply(seq_len(q), function(j) {
      vj <- v[[j]]
      h <- function(b) (b - vj) / step - ly(b)
      g <- G[, j]
      hj <- H[, j]
      cand <- unique(g[which(hj == 0)])
      for (k in which(hj[-n_grid] * hj[-1L] < 0)) {
        cand <- c(cand, .prox_uniroot(h, g[[k]], g[[k + 1L]],
                                      hj[[k]], hj[[k + 1L]]))
      }
      if (!length(cand)) return(.prox_root_widened(pen, h, vj))
      obj <- (cand - vj)^2 / (2 * step) -
        as.numeric(distributions7::distrib_pdf(pen@parent, cand, theta,
                                               log = TRUE))
      cand[[which.min(obj)]]
    }, numeric(1))
  }

  h <- heavy_penalty(n_coef = 1)
  th0 <- align_ptheta(h, list(sigma = 1, nu = 4))
  expect_identical(.prox_mode(function(b) {
    as.numeric(distributions7::distrib_grad_y(h@parent, b, th0)) + 0 * b
  }), 0)

  # the cases of the earlier repair of this route
  cs <- rbind(c(12, 8, 1, 4), c(3, 1, 0.3, 1), c(-3, 1, 0.3, 1),
              c(12, 8, 0.3, 2), c(-12, 8, 0.3, 2), c(40, 40, 0.3, 4),
              c(-40, 40, 0.3, 4), c(-12, 8, 1, 4), c(0.5, 40, 0.3, 4),
              c(-40, 0.1, 3, 30), c(12, 1, 1, 1), c(1e-11, 1, 1, 4),
              c(-1e-20, 1, 1, 4), c(1e-200, 1, 1, 4))
  for (i in seq_len(nrow(cs))) {
    th <- list(sigma = cs[i, 3], nu = cs[i, 4])
    expect_identical(penalty_prox(h, cs[i, 1], cs[i, 2], th),
                     scan0(h, cs[i, 1], cs[i, 2], th))
  }
  # a heavy_penalty vector
  set.seed(3)
  h200 <- heavy_penalty(n_coef = 200)
  x <- stats::rnorm(200, sd = 3)
  for (a in list(list(1, list(sigma = 1, nu = 4)),
                 list(7, list(sigma = 0.2, nu = 1.5)))) {
    expect_identical(penalty_prox(h200, x, a[[1]], a[[2]]),
                     scan0(h200, x, a[[1]], a[[2]]))
  }
  # a fixed() Student t written with its location at zero
  p0 <- t1_at(0)
  for (v in c(-30, -2, 0.4, 5, 30)) {
    expect_identical(penalty_prox(p0, v, 4, list(sigma = 0.01, nu = 1)),
                     scan0(p0, v, 4, list(sigma = 0.01, nu = 1)))
  }
})
