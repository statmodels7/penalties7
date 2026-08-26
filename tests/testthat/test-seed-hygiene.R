# check_penalty() draws a coefficient vector when beta is left at its default.
# Until 0.18.0 it called set.seed(7) and never put back what it found, so a call
# inside a simulation silently changed the simulation. The seed stays fixed, so
# a validator reports the same worst error on two runs; only the restore is new.

test_that("the default draw leaves the caller's stream alone", {
  set.seed(99)
  want <- runif(1)
  set.seed(99)
  invisible(check_penalty(lasso_penalty(n_coef = 3), verbose = FALSE))
  expect_identical(runif(1), want)
})

test_that("passing beta leaves it alone too, as it always did", {
  set.seed(99)
  want <- runif(1)
  set.seed(99)
  invisible(check_penalty(lasso_penalty(n_coef = 3), beta = c(1, 2, 3),
                          verbose = FALSE))
  expect_identical(runif(1), want)
})

test_that("the report does not depend on the state it was called from", {
  a <- check_penalty(lasso_penalty(n_coef = 4), verbose = FALSE)
  set.seed(1234)
  invisible(rnorm(31))
  b <- check_penalty(lasso_penalty(n_coef = 4), verbose = FALSE)
  expect_identical(a, b)
})

test_that("no seed is left behind when the caller had none", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(check_penalty(ridge_penalty(n_coef = 3), verbose = FALSE))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the draw itself is unchanged, so no reported number moves", {
  # the discriminator: the seed is still 7, so the beta the validator uses is
  # the one it always used and every statistic is the one it always reported
  set.seed(7)
  want <- round(stats::rnorm(3, sd = 1.3), 2) + 0.11
  p <- ridge_penalty(n_coef = 3)          # no kinks, so no push is taken
  expect_identical(check_penalty(p, verbose = FALSE),
                   check_penalty(p, beta = want, verbose = FALSE))
})
