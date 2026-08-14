# The three branches, each against a route that shares no code with it.

test_that("the quadratic penalty passes its battery and knows its rank", {
  P <- crossprod(diff(diag(5), differences = 2))
  pen <- quadratic_penalty(P)
  res <- check_penalty(pen, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
  expect_identical(penalty_rank(pen), 3L)
  expect_false(is_proper(pen))
  expect_true(is_proper(quadratic_penalty(diag(4))))
  # the null space of a second-difference penalty is the linear polynomials
  nb <- penalty_null_basis(pen)
  expect_identical(ncol(nb), 2L)
  expect_lt(max(abs(crossprod(P, nb))), 1e-10)
})

test_that("the separable branch passes for the Gaussian prior and the t", {
  gauss <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian3_distrib(), mu = 0),
    n_coef = 4)
  for (pen in list(gauss, heavy_penalty(n_coef = 4))) {
    th <- if (identical(pen@params, "tau")) list(tau = 1 / 1.4^2) else
      list(sigma = 1.4, nu = 6)
    res <- check_penalty(pen, theta = th, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(pen@penalty_name,
                             paste(res$check[res$status != "OK"],
                                   collapse = ", ")))
  }
  # the lasso away from its declared kink
  pen <- lasso_penalty(n_coef = 3)
  expect_identical(penalty_kinks(pen, list(lambda = 1)), 0)
  res <- check_penalty(pen, beta = c(0.7, -1.2, 2.1), theta = list(lambda = 1),
                       verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("ridge exists twice and the two constructions agree exactly", {
  # quadratic_penalty(I) at lambda = 1/sigma^2 IS the gaussian at scale
  # sigma, up to the same constant: two implementations of one object, so
  # no tolerance has to be chosen beyond machine precision
  q <- 3
  beta <- c(0.4, -1.1, 2.2)
  sigma <- 1.7
  a <- quadratic_penalty(diag(q))
  b <- distrib_penalty(
    distributions7::fixed(distributions7::gaussian3_distrib(), mu = 0),
    n_coef = q)
  tha <- list(lambda = 1 / sigma^2)
  thb <- list(tau = 1 / sigma^2)
  expect_equal(penalty_value(a, beta, tha), penalty_value(b, beta, thb),
               tolerance = 1e-14)
  expect_equal(penalty_gradient(a, beta, tha), penalty_gradient(b, beta, thb),
               tolerance = 1e-14)
  expect_equal(penalty_hessian(a, beta, tha), penalty_hessian(b, beta, thb),
               tolerance = 1e-14)
})

test_that("SCAD and MCP match the literature and pass the battery", {
  pen <- scad_penalty(n_coef = 6)
  th <- list(lambda = 1, a = 3.7)
  # rho' typed independently from Fan-Li, on a grid inside the regions
  t <- c(0.3, 0.8, 1.5, 3.0, 4.5, 8.0)
  d1 <- penalty_gradient(pen, t, th)
  ref <- ifelse(t <= 1, 1, ifelse(t <= 3.7, (3.7 - t) / 2.7, 0))
  expect_equal(d1, ref, tolerance = 1e-14)
  expect_equal(penalty_value(pen, 0 * t, th), 0)
  # constant beyond the shoulder
  expect_equal(penalty_value(pen, 10, th), penalty_value(pen, 100, th))
  res <- check_penalty(pen, beta = t, theta = th, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))

  pen <- mcp_penalty(n_coef = 4)
  th <- list(lambda = 1, gamma = 3)
  t <- c(0.4, 1.2, 2.4, 5.0)
  d1 <- penalty_gradient(pen, t, th)
  expect_equal(d1, pmax(1 - t / 3, 0) * 1, tolerance = 1e-14)
  res <- check_penalty(pen, beta = t, theta = th, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("the map enters every derivative the same way", {
  D <- diff(diag(4))
  pen <- lasso_penalty(map = D)   # the fused lasso
  th <- list(lambda = 0.8)
  beta <- c(0.5, 1.4, -0.3, 2.0)
  res <- check_penalty(pen, beta = beta, theta = th, verbose = FALSE)
  expect_true(all(res$status == "OK"))
  # gradient through the map equals D' applied to the identity-map gradient
  flat <- lasso_penalty(n_coef = 3)
  g0 <- penalty_gradient(flat, as.numeric(D %*% beta), th)
  expect_equal(penalty_gradient(pen, beta, th),
               as.numeric(crossprod(D, g0)), tolerance = 1e-14)
})

test_that("the link scale is the chain rule on the parameter scale", {
  skip_if_not_installed("numDeriv")
  pen <- quadratic_penalty(diag(3))
  beta <- c(1, -0.5, 0.2)
  th <- list(lambda = 2)
  gl <- penalty_grad_theta(pen, beta, th, scale = "link")[[1]]
  ref <- numDeriv::grad(function(e) {
    penalty_value(pen, beta, list(lambda = exp(e)))
  }, log(2))
  expect_equal(gl, ref, tolerance = 1e-8)
  hl <- penalty_hess_theta(pen, beta, th, scale = "link")[[1]]
  refh <- numDeriv::hessian(function(e) {
    penalty_value(pen, beta, list(lambda = exp(e)))
  }, log(2))
  expect_equal(hl, as.numeric(refh), tolerance = 1e-6)
})

test_that("an injected defect is caught and a clean reference still passes", {
  pen <- quadratic_penalty(diag(3))
  res <- check_penalty(pen, verbose = FALSE)
  expect_true(all(res$status == "OK"))
  # a 5% wrong gradient must fail: probe through a corrupted subclass method
  Bad <- S7::new_class("Bad", parent = QuadraticPenalty)
  S7::method(penalty_gradient, Bad) <- function(pen, beta, theta, ...) {
    1.05 * theta[[1]] * as.numeric(pen@DPD %*% beta)
  }
  bad <- do.call(Bad, S7::props(quadratic_penalty(diag(3))))
  res <- check_penalty(bad, verbose = FALSE)
  expect_true(any(res$status == "FAILED"))
})

test_that("the base class refuses the marginal pieces off the quadratic branch", {
  pen <- lasso_penalty(n_coef = 2)
  expect_false(is_quadratic(pen))
  expect_error(penalty_matrix(pen, list(lambda = 1)), "quadratic")
  expect_error(penalty_rank(pen), "quadratic")
  expect_error(penalty_logpdet(pen, list(lambda = 1)), "quadratic")

  # and the ridge IS quadratic now: a Gaussian prior written by its precision
  # is the quadratic penalty at the identity, the same value to the last bit
  r <- ridge_penalty(n_coef = 2)
  expect_true(is_quadratic(r))
  expect_identical(r@params, "lambda")
})


test_that("the structured prior passes its battery on three structures", {
  skip_if_not_installed("numDeriv")
  for (s in list(parameters7::log_cholesky(3, role = "precision"),
                 parameters7::ar1(4, role = "precision"),
                 parameters7::compound_symmetry(3, role = "precision"))) {
    pen <- structured_penalty(s)
    set.seed(5)
    th <- stats::setNames(as.list(round(rnorm(s@n_free, sd = 0.3), 2)),
                          s@free_names)
    res <- check_penalty(pen, theta = th, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(pen@penalty_name,
                             paste(res$check[res$status != "OK"],
                                   collapse = ", ")))
  }
})


test_that("the structured prior at the identity precision is the plain ridge", {
  # log_cholesky at a zero free vector is the identity matrix, so the prior
  # is quadratic_penalty(diag(q)) at lambda = 1: two constructions of one
  # object, compared with no tolerance to choose
  q <- 3
  beta <- c(0.4, -1.1, 2.2)
  s <- parameters7::log_cholesky(q, role = "precision")
  pen <- structured_penalty(s)
  th <- stats::setNames(as.list(rep(0, s@n_free)), s@free_names)
  ref <- quadratic_penalty(diag(q))
  expect_equal(penalty_value(pen, beta, th),
               penalty_value(ref, beta, list(lambda = 1)), tolerance = 1e-14)
  expect_equal(penalty_gradient(pen, beta, th),
               penalty_gradient(ref, beta, list(lambda = 1)),
               tolerance = 1e-14)
  expect_equal(unname(penalty_hessian(pen, beta, th)),
               penalty_hessian(ref, beta, list(lambda = 1)),
               tolerance = 1e-14)
  expect_true(is_quadratic(pen))
  expect_true(is_proper(pen))
  expect_identical(penalty_rank(pen), 3L)
  expect_equal(penalty_logpdet(pen, th)$value, 0, tolerance = 1e-14)
})


test_that("a named vector of hyperparameters is read like the list", {
  # The branches split on how they read theta: `[[` accepts a named numeric
  # vector as well as a list, `$` accepts only the list. A caller passing a
  # vector therefore reached the quadratic and separable branches and stopped
  # inside scad and mcp with "$ operator is invalid for atomic vectors",
  # three frames down and naming neither the argument nor the penalty.
  beta <- c(0.4, -1.1, 2.2, 0.05)
  cases <- list(
    quadratic = list(pen = quadratic_penalty(diag(4)), th = c(lambda = 2.5)),
    ridge     = list(pen = ridge_penalty(n_coef = 4L), th = c(lambda = 1.5)),
    lasso     = list(pen = lasso_penalty(n_coef = 4L), th = c(lambda = 1.7)),
    scad      = list(pen = scad_penalty(n_coef = 4L),
                     th = c(lambda = 1.2, a = 3.7)),
    mcp       = list(pen = mcp_penalty(n_coef = 4L),
                     th = c(lambda = 1.2, gamma = 3))
  )
  for (nm in names(cases)) {
    pen <- cases[[nm]]$pen
    v <- cases[[nm]]$th
    l <- as.list(v)
    expect_equal(penalty_value(pen, beta, v), penalty_value(pen, beta, l),
                 tolerance = 1e-14, label = nm)
    expect_equal(penalty_gradient(pen, beta, v),
                 penalty_gradient(pen, beta, l), tolerance = 1e-14,
                 label = nm)
    expect_equal(penalty_kinks(pen, v), penalty_kinks(pen, l), label = nm)
    expect_equal(penalty_prox(pen, beta, 0.3, v),
                 penalty_prox(pen, beta, 0.3, l), tolerance = 1e-14,
                 label = nm)
  }

  # and the reordering the alignment already did is unaffected by the shape
  p <- scad_penalty(n_coef = 4L)
  expect_equal(penalty_value(p, beta, c(a = 3.7, lambda = 1.2)),
               penalty_value(p, beta, list(lambda = 1.2, a = 3.7)),
               tolerance = 1e-14)
  # a missing hyperparameter is still reported, whichever shape it arrives in
  expect_error(penalty_value(p, beta, c(lambda = 1.2)), "Missing parameter")
})


test_that("a separable penalty derives its kinks from the parent", {
  # kinks was a constructor argument with no default beyond "none", so a
  # penalty built by hand from a non-smooth parent declared itself smooth and
  # the model layer put it in the scheme for the opposite property. The
  # candidates come from params_smooth crossed with what fixed() holds, and
  # each one is then measured: inferring alone would put a kink on any family
  # whose non-smooth parameter is not a location.
  lap <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
  expect_equal(distrib_kinks(lap), 0)
  expect_equal(distrib_kinks(distributions7::fixed(
    distributions7::laplace_distrib(), mu = 1.5)), 1.5)
  # the rate parametrization is the same law and answers the same
  expect_equal(distrib_kinks(distributions7::fixed(
    distributions7::laplace2_distrib(), mu = 0)), 0)
  # the elastic net's density carries an absolute value too
  expect_equal(distrib_kinks(distributions7::fixed(
    distributions7::enet_distrib(), mu = 0)), 0)

  # a smooth family has none, whatever is fixed
  expect_length(distrib_kinks(distributions7::fixed(
    distributions7::gaussian1_distrib(), mu = 0)), 0L)
  expect_length(distrib_kinks(distributions7::fixed(
    distributions7::student_t1_distrib(), mu = 0)), 0L)
  # and nothing is taken from a parameter that is free: where it sits is
  # whatever the hyperparameters say at the time
  expect_length(distrib_kinks(distributions7::laplace_distrib()), 0L)

  # the penalty carries what the parent implies
  pen <- distrib_penalty(lap, n_coef = 5L)
  expect_equal(penalty_kinks(pen, list(sigma = 1)), 0)
  # and a caller who says otherwise is obeyed
  expect_length(penalty_kinks(distrib_penalty(lap, n_coef = 5L,
                                              kinks = numeric(0)),
                              list(sigma = 1)), 0L)
  expect_equal(penalty_kinks(distrib_penalty(lap, n_coef = 5L, kinks = c(0, 2)),
                             list(sigma = 1)), c(0, 2))

  # the shipped instances are unchanged
  expect_equal(penalty_kinks(lasso_penalty(n_coef = 3L), list(lambda = 1)), 0)
  expect_equal(penalty_kinks(elasticnet_penalty(n_coef = 3L),
                             list(lambda = 1, alpha = 0.5)), 0)
  expect_length(penalty_kinks(ridge_penalty(n_coef = 3L), list(lambda = 1)), 0L)
})


test_that("a Matrix map does not decide the class of what a penalty returns", {
  # 0.10.0 let a map stay a Matrix, which is what makes a diagonal one cost q
  # numbers instead of q^2. The congruence D' H D then carried that class into
  # the RESULT, so a consumer writing a penalty's Hessian into a block of its
  # own information failed on the class rather than on the arithmetic -- which
  # is how a standardized term first broke a fit. Every branch, every route.
  s <- c(0.5, 2, 4)
  b <- c(0.3, -1.2, 0.7)
  P <- crossprod(diff(diag(3)))
  maps <- list(identity = NULL, diagonal = Matrix::Diagonal(x = s),
               dense = diag(s))
  for (tag in names(maps)) {
    map <- maps[[tag]]
    pens <- list(
      ridge = ridge_penalty(map = map, n_coef = 3L),
      lasso = lasso_penalty(map = map, n_coef = 3L),
      enet  = elasticnet_penalty(map = map, n_coef = 3L),
      scad  = scad_penalty(map = map, n_coef = 3L),
      mcp   = mcp_penalty(map = map, n_coef = 3L),
      quadratic = quadratic_penalty(P, map = map),
      additive  = additive_penalty(list(P, diag(3)), map = map)
    )
    for (nm in names(pens)) {
      pen <- pens[[nm]]
      th <- lapply(pen@params_bounds, function(bb) {
        if (all(is.finite(bb))) mean(bb) else if (is.finite(bb[1])) bb[1] + 1
        else 1
      })
      lab <- paste(tag, nm)
      expect_true(is.numeric(penalty_value(pen, b, th)), label = lab)
      expect_true(is.numeric(penalty_gradient(pen, b, th)), label = lab)
      h <- penalty_hessian(pen, b, th)
      expect_true(is.matrix(h) && is.numeric(h), label = lab)
      # and a base matrix takes the subassignment every consumer writes
      m <- matrix(0, 3, 3)
      expect_silent(m[1:3, 1:3] <- m[1:3, 1:3] + h)
    }
  }
})


# ---------------------------------------------------------------------------
# the structured prior reads its structure's role
# ---------------------------------------------------------------------------

test_that("a structure that does not say which matrix it is is rejected", {
  # "either" is a statement about the structure, not about this prior, and the
  # two readings differ in the sign of the log-determinant term: a guess would
  # give a fit converging to a different matrix without saying so
  expect_error(structured_penalty(parameters7::log_cholesky(2)), "role")
  expect_error(structured_penalty(parameters7::ar1(3)), "role")
  expect_silent(structured_penalty(parameters7::ar1(3, role = "covariance")))
})


test_that("a rank-deficient covariance is rejected and a precision is not", {
  d <- parameters7::scaled_matrix(diag(c(1, 1, 0)), role = "covariance")
  expect_error(structured_penalty(d), "no inverse")
  ok <- parameters7::scaled_matrix(diag(c(1, 1, 0)), role = "precision")
  expect_false(is_proper(structured_penalty(ok)))
})


test_that("the covariance branch at Sigma is the precision branch at Sigma^-1", {
  # the twin that makes the SIGN verifiable without choosing a tolerance: the
  # two are the same prior, so every quantity agrees to the last bit
  set.seed(11)
  for (d in 2:4) {
    sc <- parameters7::log_cholesky(d, role = "covariance")
    sp <- parameters7::log_cholesky(d, role = "precision")
    eta <- stats::rnorm(sc@n_free, sd = 0.4)
    b <- stats::rnorm(d)
    om <- solve(unclass(parameters7::param_value(sc, eta)))
    th_c <- stats::setNames(as.list(eta), sc@free_names)
    th_p <- stats::setNames(as.list(parameters7::param_free(sp, om)),
                            sp@free_names)
    pc <- structured_penalty(sc)
    pp <- structured_penalty(sp)
    expect_equal(penalty_value(pc, b, th_c), penalty_value(pp, b, th_p),
                 tolerance = 1e-12, info = paste("d =", d))
    expect_equal(penalty_gradient(pc, b, th_c), penalty_gradient(pp, b, th_p),
                 tolerance = 1e-12)
    expect_equal(unname(penalty_hessian(pc, b, th_c)),
                 unname(penalty_hessian(pp, b, th_p)), tolerance = 1e-12)
    # and the value IS the negative log-density of N(0, Sigma)
    sig <- unclass(parameters7::param_value(sc, eta))
    ref <- -(-d / 2 * log(2 * pi) -
             determinant(sig, logarithm = TRUE)$modulus[1] / 2 -
             sum(b * solve(sig, b)) / 2)
    expect_equal(penalty_value(pc, b, th_c), ref, tolerance = 1e-12)
  }
})


test_that("the covariance branch's theta derivatives agree with numDeriv", {
  skip_if_not_installed("numDeriv")
  # d = 2 says nothing about the cross terms, so the loop starts where the
  # second derivative in two DIFFERENT free values exists
  set.seed(12)
  for (d in 2:4) {
    s <- parameters7::log_cholesky(d, role = "covariance")
    pen <- structured_penalty(s)
    eta <- stats::rnorm(s@n_free, sd = 0.35)
    b <- stats::rnorm(d)
    nm <- s@free_names
    fv <- function(z) penalty_value(pen, b, stats::setNames(as.list(z), nm))
    th <- stats::setNames(as.list(eta), nm)
    expect_equal(unname(unlist(penalty_grad_theta(pen, b, th))),
                 numDeriv::grad(fv, eta), tolerance = 1e-6,
                 info = paste("d =", d))
    hh <- penalty_hess_theta(pen, b, th)
    H <- numDeriv::hessian(fv, eta)
    for (i in seq_along(nm)) for (j in seq_along(nm)) {
      key <- paste(nm[sort(c(i, j))], collapse = "_")
      if (i == j) key <- paste0(nm[i], "_", nm[i])
      if (!is.null(hh[[key]])) {
        expect_equal(hh[[key]], H[i, j], tolerance = 1e-5,
                     info = paste("d =", d, key))
      }
    }
    fg <- function(z) penalty_gradient(pen, b, stats::setNames(as.list(z), nm))
    expect_equal(unname(do.call(cbind, penalty_cross(pen, b, th))),
                 numDeriv::jacobian(fg, eta), tolerance = 1e-6)
  }
})


test_that("the covariance branch passes check_penalty", {
  skip_if_not_installed("numDeriv")
  for (s in list(parameters7::log_cholesky(3, role = "covariance"),
                 parameters7::ar1(4, role = "covariance"),
                 parameters7::compound_symmetry(3, role = "covariance"))) {
    pen <- structured_penalty(s)
    set.seed(5)
    th <- stats::setNames(as.list(round(stats::rnorm(s@n_free, sd = 0.3), 2)),
                          s@free_names)
    res <- check_penalty(pen, theta = th, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(pen@penalty_name,
                             paste(res$check[res$status != "OK"],
                                   collapse = ", ")))
  }
})
