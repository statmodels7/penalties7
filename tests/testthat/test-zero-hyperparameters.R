# A penalty with no free hyperparameters is a legitimate object:
# distributions7::fixed() documents n_params = 0 as legal, and it is what a
# caller builds to hold a prior at a value the outer search must not touch.
# Until 0.17.0 three generics stopped on it, two of them through a paste0
# recycling in ptheta_pairs() and one through which() of an empty comparison.

d7 <- asNamespace("distributions7")

held_gaussian <- function() {
  distrib_penalty(d7$fixed(d7$gaussian1_distrib(), mu = 0, sigma = 1), n_coef = 3)
}

test_that("ptheta_pairs() answers the empty named list for no hyperparameters", {
  # named, so penalty_hess_theta() keeps the shape its two siblings already
  # had for this case
  expect_identical(ptheta_pairs(character(0)), stats::setNames(list(), character(0)))
  expect_identical(names(ptheta_pairs(character(0))), character(0))
  # and still keys one and two the way it did
  expect_identical(names(ptheta_pairs("lambda")), "lambda_lambda")
  expect_identical(names(ptheta_pairs(c("lambda", "alpha"))),
                   c("lambda_lambda", "alpha_alpha", "lambda_alpha"))
})

test_that("the guard is what fixes it: paste0 still recycles without it", {
  # the negative control
  expect_identical(paste0(character(0), "_", character(0)), "_")
  expect_error(stats::setNames(list(), "_"))
})

test_that("every generic answers for a penalty with no hyperparameters", {
  p <- held_gaussian()
  b <- c(0.4, -1.1, 0.7)
  th <- list()
  expect_length(p@params, 0L)

  expect_type(penalty_value(p, b, th), "double")
  expect_length(penalty_gradient(p, b, th), 3L)
  expect_true(is.matrix(penalty_hessian(p, b, th)))
  # all three theta blocks answer, and all three answer with the same shape
  empty <- stats::setNames(list(), character(0))
  expect_identical(penalty_grad_theta(p, b, th), empty)
  expect_identical(penalty_hess_theta(p, b, th), empty)
  expect_identical(penalty_cross(p, b, th), empty)
  expect_length(penalty_prox(p, b, 0.5, th), 3L)

  r <- check_penalty(p, beta = b, verbose = FALSE)
  expect_false(any(r$status == "FAIL"))
})

test_that("holding a parameter gives the same answer as leaving it free", {
  # the discriminating comparison: a held scale is the same penalty as a free
  # one read at that scale, and every generic already agreed except the prox
  b <- c(0.4, -1.1, 0.7)
  step <- 0.5
  cases <- list(
    list(held = d7$fixed(d7$gaussian1_distrib(), mu = 0, sigma = 1),
         free = d7$fixed(d7$gaussian1_distrib(), mu = 0),
         theta = list(sigma = 1)),
    list(held = d7$fixed(d7$laplace2_distrib(), mu = 0, lambda = 1.3),
         free = d7$fixed(d7$laplace2_distrib(), mu = 0),
         theta = list(lambda = 1.3)),
    list(held = d7$fixed(d7$enet_distrib(), mu = 0, lambda = 2, alpha = 0.6),
         free = d7$fixed(d7$enet_distrib(), mu = 0),
         theta = list(lambda = 2, alpha = 0.6))
  )
  for (cs in cases) {
    ph <- distrib_penalty(cs$held, n_coef = 3)
    pf <- distrib_penalty(cs$free, n_coef = 3)
    expect_equal(penalty_prox(ph, b, step, list()),
                 penalty_prox(pf, b, step, cs$theta), tolerance = 0)
    expect_equal(penalty_value(ph, b, list()),
                 penalty_value(pf, b, cs$theta), tolerance = 0)
    # the table the compiled coordinate descent reads takes the same route
    expect_equal(penalty_prox_spec(ph, list(), step),
                 penalty_prox_spec(pf, cs$theta, step), tolerance = 0)
  }
})

test_that("one held and one free is the same answer again", {
  b <- c(0.4, -1.1, 0.7)
  pm <- distrib_penalty(d7$fixed(d7$enet_distrib(), mu = 0, lambda = 2), n_coef = 3)
  pf <- distrib_penalty(d7$fixed(d7$enet_distrib(), mu = 0), n_coef = 3)
  expect_identical(pm@params, "alpha")
  expect_equal(penalty_prox(pm, b, 0.5, list(alpha = 0.6)),
               penalty_prox(pf, b, 0.5, list(lambda = 2, alpha = 0.6)),
               tolerance = 0)
})

test_that(".prox_param names a parameter the parent does not carry", {
  p <- held_gaussian()
  expect_error(.prox_param(p, list(), "nonesuch"), "no parameter 'nonesuch'")
})
