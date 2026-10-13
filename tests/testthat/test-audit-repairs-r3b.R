# The interpretable quantities of a structured penalty, read from the
# structure through parameters7::param_readable().

test_that("an ar1 structured penalty reports the scale and rho of its matrix", {
  pen <- structured_penalty(parameters7::ar1(4))
  th <- list(log_scale = 0.2, z_rho = 0.5)
  rd <- penalty_readable(pen, th)
  expect_false(is.null(rd))
  expect_identical(names(rd$value), c("precision_scale", "precision_rho"))

  # the matrix built here from the structure's own value: the scale is its
  # diagonal entry and the correlation the ratio of an entry next to the
  # diagonal to the diagonal entry, at every position of the band
  M <- unclass(parameters7::param_value(parameters7::ar1(4), c(0.2, 0.5)))
  sc <- unname(rd$value[["precision_scale"]])
  rho <- unname(rd$value[["precision_rho"]])
  expect_equal(sc, M[1, 1], tolerance = 1e-14)
  expect_equal(sc, M[4, 4], tolerance = 1e-14)
  for (i in 1:3) {
    expect_equal(rho, M[i, i + 1] / M[i, i], tolerance = 1e-13)
  }
  # and at lag two the square, which is what makes it the AR(1) correlation
  expect_equal(rho^2, M[1, 3] / M[1, 1], tolerance = 1e-13)

  # the matrix is the precision of the prior, the one the penalty uses
  expect_equal(unclass(penalty_matrix(pen, th)), M, tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_identical(unname(rd$transform), c("log", "atanh"))
})


test_that("a compound symmetry structured penalty reports its scale and rho", {
  pen <- structured_penalty(parameters7::compound_symmetry(3))
  th <- stats::setNames(list(-0.3, 0.4), pen@params)
  rd <- penalty_readable(pen, th)
  expect_false(is.null(rd))
  expect_identical(names(rd$value), c("precision_scale", "precision_rho"))

  M <- unclass(parameters7::param_value(parameters7::compound_symmetry(3),
                                        c(-0.3, 0.4)))
  expect_equal(unname(rd$value[["precision_scale"]]), M[2, 2],
               tolerance = 1e-14)
  off <- M[upper.tri(M)] / M[1, 1]
  expect_equal(rep(unname(rd$value[["precision_rho"]]), 3), off,
               tolerance = 1e-13)
})


test_that("the Jacobian of a structured penalty's quantities is numDeriv's", {
  skip_if_not_installed("numDeriv")
  cases <- list(
    list(s = parameters7::ar1(4), v = c(0.2, 0.5)),
    list(s = parameters7::compound_symmetry(3), v = c(-0.3, 0.4)),
    list(s = parameters7::autoregressive(5, 2), v = c(0.1, 0.9, -0.4))
  )
  for (cs in cases) {
    pen <- structured_penalty(cs$s)
    th <- stats::setNames(as.list(cs$v), pen@params)
    rd <- penalty_readable(pen, th)
    expect_false(is.null(rd))
    f <- function(v) {
      penalty_readable(pen, stats::setNames(as.list(v), pen@params))$value
    }
    num <- numDeriv::jacobian(f, cs$v)
    expect_identical(dim(rd$jacobian), dim(num))
    expect_equal(unname(rd$jacobian), num, tolerance = 1e-8)
  }
})


test_that("a structured penalty's quantities have the separable branch's shape", {
  mv <- distrib_penalty(
    distributions7::fixed(distributions7::mvgaussian1_distrib(2),
                          mu1 = 0, mu2 = 0), n_coef = 6)
  ref <- penalty_readable(mv, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
                                   sigma_L2.1 = 0.5))

  pen <- structured_penalty(parameters7::ar1(4))
  th <- list(log_scale = 0.2, z_rho = 0.5)
  rd <- penalty_readable(pen, th)
  expect_identical(names(rd), names(ref))
  for (el in names(ref)) {
    expect_identical(typeof(rd[[el]]), typeof(ref[[el]]), info = el)
  }
  expect_true(is.matrix(rd$jacobian))
  expect_identical(dimnames(rd$jacobian), list(names(rd$value), pen@params))
  expect_identical(names(rd$transform), names(rd$value))
  expect_identical(names(rd$block), names(rd$value))
  expect_length(rd$block, length(rd$value))

  # every name says that the quantity belongs to the precision
  for (s in list(parameters7::ar1(4), parameters7::compound_symmetry(3),
                 parameters7::autoregressive(5, 2),
                 parameters7::scaled_matrix(diag(3)))) {
    sp <- structured_penalty(s)
    r <- penalty_readable(sp, stats::setNames(as.list(rep(0.1, s@n_free)),
                                              sp@params))
    expect_true(all(startsWith(names(r$value), "precision_")))
    expect_identical(rownames(r$jacobian), names(r$value))
    expect_identical(names(r$transform), names(r$value))
    expect_identical(names(r$block), names(r$value))
  }

  # theta as a named numeric vector gives the same answer
  expect_identical(penalty_readable(pen, unlist(th)), rd)

  # a structure that declares nothing gives NULL, on either side
  lc <- structured_penalty(parameters7::log_cholesky(2))
  expect_null(penalty_readable(lc, list(log_L1 = 0.1, log_L2 = -0.1,
                                        L2.1 = 0.3)))
  inv <- structured_penalty(parameters7::inverse_of(parameters7::ar1(4)))
  expect_null(penalty_readable(inv, th))
})
