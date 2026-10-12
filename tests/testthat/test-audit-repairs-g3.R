fixed0 <- function(d) distributions7::fixed(d, mu = 0)

grid_mass <- function(pen, theta, lim, h) {
  b <- seq(-lim, lim, by = h)
  s <- 0
  for (x in b) for (y in b) s <- s + exp(-penalty_value(pen, c(x, y), theta))
  s * h^2
}


test_that("is_proper() of a separable penalty is TRUE only for a square nonsingular map", {
  g <- fixed0(distributions7::gaussian1_distrib())
  expect_true(is_proper(distrib_penalty(g, n_coef = 2)))
  expect_true(is_proper(distrib_penalty(g, map = matrix(c(2, 1, 0, 3), 2))))
  expect_true(is_proper(distrib_penalty(g, map = Matrix::Diagonal(x = c(2, 3)))))
  # rank one: flat along (1, -1)
  expect_false(is_proper(distrib_penalty(g, map = matrix(1, 2, 2))))
  # tall and of full rank: the missing constant depends on sigma
  expect_false(is_proper(distrib_penalty(g, map = matrix(c(1, 0, 1, 0, 1, 1), 3))))
  # a diagonal with a zero leaves a coordinate unpenalized
  expect_false(is_proper(distrib_penalty(
    g, map = Matrix::Diagonal(x = c(1, 0, 2)))))
})


test_that("under a square nonsingular map exp(-value) integrates to one", {
  D <- matrix(c(2, 1, 0, 3), 2)  # det 6, not diagonal
  g <- distrib_penalty(fixed0(distributions7::gaussian1_distrib()), map = D)
  expect_equal(grid_mass(g, list(sigma = 1), lim = 4, h = 0.1), 1,
               tolerance = 1e-8)
  # the kinks of the Laplace sit off the grid, so the trapezoid error is
  # O(h^2): 4e-3 at this width
  l <- distrib_penalty(fixed0(distributions7::laplace2_distrib()), map = D)
  expect_equal(grid_mass(l, list(lambda = 2), lim = 5, h = 0.05), 1,
               tolerance = 0.02)
  # a bivariate parent read as one block: the Jacobian is the whole map's
  mv <- distrib_penalty(
    distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                          mu1 = 0, mu2 = 0), map = D)
  expect_equal(grid_mass(mv, list(sigma_log_L1 = 0, sigma_log_L2 = 0,
                                  sigma_L2.1 = 0), lim = 4, h = 0.1), 1,
               tolerance = 1e-8)

  expect_equal(g@log_jacobian, log(6), tolerance = 1e-14)

  # a diagonal map keeps the value of the earlier formula to the last bit,
  # and a diagonal with a zero still subtracts nothing
  b <- c(0.7, -1.3, 0.2)
  for (dg in list(c(2, 0.5, 3), c(1, 0, 2))) {
    for (pp in list(list(fixed0(distributions7::laplace2_distrib()),
                         list(lambda = 1.7)),
                    list(fixed0(distributions7::student_t1_distrib()),
                         list(sigma = 0.8, nu = 3)))) {
      pen <- distrib_penalty(pp[[1]], map = Matrix::Diagonal(x = dg))
      v <- -sum(distributions7::distrib_pdf(pp[[1]], dg * b, pp[[2]],
                                            log = TRUE))
      if (all(dg != 0)) v <- v - sum(log(abs(dg)))
      expect_identical(penalty_value(pen, b, pp[[2]]), v)
    }
  }
})


test_that("a diagonal base matrix is the same map as a Matrix diagonal", {
  par <- fixed0(distributions7::laplace2_distrib())
  pa <- distrib_penalty(par, map = diag(c(2, 3)))
  pb <- distrib_penalty(par, map = Matrix::Diagonal(x = c(2, 3)))
  b <- c(0.4, -1.1)
  th <- list(lambda = 1.3)
  expect_s4_class(pa@map, "diagonalMatrix")
  expect_identical(penalty_value(pa, b, th), penalty_value(pb, b, th))
  expect_identical(has_prox(pa), has_prox(pb))
  expect_true(has_prox(pa))
  expect_identical(penalty_prox(pa, b, 0.5, th), penalty_prox(pb, b, 0.5, th))
  set.seed(11); da <- penalty_draw(pa, th)
  set.seed(11); db <- penalty_draw(pb, th)
  expect_identical(da, db)
  expect_true(all(is.finite(da)))
  # a dense Matrix that is diagonal is converted too, and an off-diagonal
  # entry keeps the map as it is
  pc <- distrib_penalty(par, map = Matrix::Matrix(diag(c(2, 3)), sparse = FALSE))
  expect_s4_class(pc@map, "diagonalMatrix")
  expect_false(methods::is(
    distrib_penalty(par, map = matrix(c(2, 1e-300, 0, 3), 2))@map,
    "diagonalMatrix"))
})


test_that("a discrete parent is rejected, naming it", {
  bad <- list(distributions7::poisson_distrib(),
              distributions7::bernoulli_distrib(),
              distributions7::fixed(distributions7::poisson_distrib(), mu = 1),
              distributions7::zero_inflated(distributions7::poisson_distrib()),
              distributions7::multinomial_distrib(n_dim = 2, size = 5))
  for (d in bad) {
    expect_error(distrib_penalty(d, n_coef = 2), d@distrib_name, fixed = TRUE)
    expect_error(distrib_penalty(d, n_coef = 2), "continuous density")
  }
  ok <- list(fixed0(distributions7::gaussian1_distrib()),
             distributions7::truncated(fixed0(distributions7::gaussian1_distrib()),
                                       lower = -1),
             distributions7::mvgaussian1_distrib(n_dim = 2))
  for (d in ok) {
    expect_true(S7::S7_inherits(distrib_penalty(d, n_coef = 2), DistribPenalty))
  }
})


test_that("has_jump() does not read a smooth point as a kink", {
  g <- fixed0(distributions7::gaussian1_distrib())
  for (s in c(0.1, 0.5, 1, 2)) {
    expect_false(has_jump(g, list(sigma = s), 0))
    expect_false(has_jump(g, list(sigma = s), 1.3))
  }
  l2 <- fixed0(distributions7::laplace2_distrib())
  expect_true(has_jump(l2, list(lambda = 1), 0))
  en <- fixed0(distributions7::enet_distrib())
  expect_true(has_jump(en, probe_theta(en), 0))
  lp <- fixed0(distributions7::laplace_distrib())
  expect_false(has_jump(lp, list(sigma = 1), 1.3))
})


test_that("distrib_kinks() finds a kink under a wrapper around fixed()", {
  l2 <- fixed0(distributions7::laplace2_distrib())
  tr <- distributions7::truncated(l2, lower = -1, upper = 3)
  # the slope of the log-density jumps by 2 lambda = 2 at zero
  th <- probe_theta(tr)
  jump <- as.numeric(distributions7::distrib_grad_y(tr, 1e-3, th)) -
    as.numeric(distributions7::distrib_grad_y(tr, -1e-3, th))
  expect_equal(jump, -2, tolerance = 1e-10)
  expect_equal(distrib_kinks(tr), 0)
  pen <- distrib_penalty(tr, n_coef = 3)
  expect_equal(pen@kinks, 0)
  expect_equal(penalty_kinks(pen, list(lambda = 1)), 0)
  # the reverse nesting, which was already found
  expect_equal(distrib_kinks(distributions7::fixed(
    distributions7::truncated(distributions7::laplace2_distrib(),
                              lower = -1, upper = 3), mu = 0)), 0)
  # a smooth parent under the same wrapper has none
  expect_length(distrib_kinks(distributions7::truncated(
    fixed0(distributions7::gaussian1_distrib()), lower = -1, upper = 3)), 0L)
  # a wrapper that moves the point drops the candidate at the confirmation
  af <- distributions7::transformation(
    l2, distributions7::affine_transform(loc = 1.5, scale = 2))
  expect_length(distrib_kinks(af), 0L)
})


test_that("params_smooth of a separable penalty is the parent's", {
  expect_identical(
    distrib_penalty(distributions7::laplace_distrib(), n_coef = 3)@params_smooth,
    c(mu = FALSE, sigma = TRUE))
  expect_identical(
    distrib_penalty(fixed0(distributions7::laplace2_distrib()),
                    n_coef = 3)@params_smooth,
    c(lambda = TRUE))
})


test_that("penalty_draw() under a diagonal with a zero is NA only there", {
  pen <- lasso_penalty(map = Matrix::Diagonal(x = c(1, 0, 2)))
  th <- list(lambda = 1.5)
  set.seed(5)
  out <- penalty_draw(pen, th)
  set.seed(5)
  t <- as.numeric(distributions7::distrib_rng(pen@parent, 3L, th))
  expect_true(is.na(out[2]))
  expect_true(all(is.finite(out[c(1, 3)])))
  expect_identical(out[c(1, 3)], t[c(1, 3)] / c(1, 2))
  # a map that is not diagonal still gives no draw
  expect_true(all(is.na(penalty_draw(
    lasso_penalty(map = matrix(c(2, 1, 0, 3), 2)), th))))
})


test_that("has_jump() finds the kink of a family written with its density only", {
  # The slope of a density-only family is differenced by the parent's own
  # fallback, which smears the jump over a few of its steps; the kink must
  # still be resolved, and a smooth density-only family must still be read as
  # smooth, far from the origin included.
  pdf <- distributions7::distrib_pdf
  UserLap <- S7::new_class("UserLapP7", parent = distributions7::continuous_distrib,
                           package = NULL)
  S7::method(pdf, UserLap) <- function(distrib, y, theta, log = FALSE, ...) {
    v <- -abs(y - theta[[1]]) / theta[[2]] - log(2 * theta[[2]])
    if (log) v else exp(v)
  }
  UserG <- S7::new_class("UserGP7", parent = distributions7::continuous_distrib,
                         package = NULL)
  S7::method(pdf, UserG) <- function(distrib, y, theta, log = FALSE, ...) {
    stats::dnorm(y, theta[[1]], theta[[2]], log = log)
  }
  make <- function(cls, nm) cls(
    distrib_name = nm, dimension = "univariate", bounds = c(-Inf, Inf),
    params = c("mu", "sigma"),
    params_interpretation = c(mu = "location", sigma = "scale"), n_params = 2,
    params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
    link_params = list(mu = linkfunctions7::identity_link(),
                       sigma = linkfunctions7::log_link()),
    params_smooth = c(mu = FALSE, sigma = TRUE))
  ul <- make(UserLap, "user laplace")
  ug <- make(UserG, "user gauss")
  for (m in c(0, 1.5, 20)) {
    expect_equal(distrib_kinks(distributions7::fixed(ul, mu = m)), m)
    expect_length(distrib_kinks(distributions7::fixed(ug, mu = m)), 0L)
  }
  g <- fixed0(distributions7::gaussian1_distrib())
  expect_false(has_jump(g, list(sigma = 0.1), 20))
})


test_that("the edge of a folded parent's support is not a kink", {
  l2 <- fixed0(distributions7::laplace2_distrib())
  fo <- distributions7::folded(l2)
  expect_length(distrib_kinks(fo), 0L)
  expect_length(distrib_penalty(fo, n_coef = 2)@kinks, 0L)
  # an interior point of a truncated support is still one
  expect_equal(distrib_kinks(distributions7::truncated(l2, lower = -1,
                                                       upper = 3)), 0)
})


test_that("penalty_null_basis() of a separable penalty is the null space of its map", {
  g <- fixed0(distributions7::gaussian1_distrib())
  expect_identical(dim(penalty_null_basis(distrib_penalty(g, n_coef = 3))),
                   c(3L, 0L))
  # a tall map of full column rank leaves no direction free
  tall <- distrib_penalty(g, map = matrix(c(1, 0, 1, 0, 1, 1), 3))
  expect_identical(dim(penalty_null_basis(tall)), c(2L, 0L))
  # a zero on a diagonal map frees exactly that coordinate
  z <- penalty_null_basis(lasso_penalty(map = Matrix::Diagonal(x = c(1, 0, 2))))
  expect_equal(z, matrix(c(0, 1, 0), 3, 1))
  # a short map leaves the directions its rows do not see
  D <- diff(diag(4), differences = 2)
  nb <- penalty_null_basis(distrib_penalty(g, map = D))
  expect_identical(dim(nb), c(4L, 2L))
  expect_lt(max(abs(D %*% nb)), 1e-12)
  expect_equal(crossprod(nb), diag(2), tolerance = 1e-12)
  # a square map of rank one leaves one
  nb1 <- penalty_null_basis(distrib_penalty(g, map = matrix(1, 2, 2)))
  expect_equal(abs(as.numeric(nb1)), rep(1 / sqrt(2), 2), tolerance = 1e-12)
  # the value does not move along the null space
  b <- c(0.3, -1.2, 0.8, 2.1)
  pen <- distrib_penalty(g, map = D)
  expect_equal(penalty_value(pen, b + 5 * nb[, 1], list(sigma = 1)),
               penalty_value(pen, b, list(sigma = 1)), tolerance = 1e-12)
})
