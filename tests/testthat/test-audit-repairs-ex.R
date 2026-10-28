test_that("SCAD and MCP store a diagonal base matrix as a diagonal map", {
  d <- c(1, 2, 3)
  b <- c(0.3, -1.2, 0.7)
  v <- c(1.5, -3, 6)
  t <- 0.2
  cases <- list(
    scad = list(ctor = scad_penalty, th = list(lambda = 1, a = 3.7)),
    mcp  = list(ctor = mcp_penalty, th = list(lambda = 1, gamma = 3))
  )
  for (nm in names(cases)) {
    for (cv in list(NULL, c(1, 1.5, 0.5))) {
      lab <- sprintf("%s, curvature %s", nm, if (is.null(cv)) "NULL" else "given")
      ctor <- cases[[nm]]$ctor
      th <- cases[[nm]]$th
      base <- ctor(map = diag(d), curvature = cv)
      mat <- ctor(map = Matrix::Diagonal(x = d), curvature = cv)
      expect_true(methods::is(base@map, "diagonalMatrix"), label = lab)
      expect_true(has_prox(base), label = lab)
      expect_identical(has_prox(base), has_prox(mat), label = lab)
      expect_equal(penalty_value(base, b, th), penalty_value(mat, b, th),
                   tolerance = 0, label = lab)
      expect_equal(penalty_gradient(base, b, th),
                   penalty_gradient(mat, b, th), tolerance = 0, label = lab)
      expect_equal(penalty_prox(base, v, t, th), penalty_prox(mat, v, t, th),
                   tolerance = 0, label = lab)
      sp <- penalty_prox_spec(base, th, t)
      expect_false(is.null(sp), label = lab)
      expect_equal(sp, penalty_prox_spec(mat, th, t), label = lab)
    }
  }
})

test_that("penalty_prox() rejects a step that is not finite and positive", {
  v <- c(0.5, -2)
  pens <- list(lasso = list(pen = lasso_penalty(n_coef = 2),
                            th = list(lambda = 1)),
               ridge = list(pen = ridge_penalty(n_coef = 2),
                            th = list(lambda = 1)))
  for (nm in names(pens)) {
    for (bad in list(Inf, NA, -1, 0, c(1, 2))) {
      expect_error(penalty_prox(pens[[nm]]$pen, v, bad, pens[[nm]]$th),
                   "'step' must be a single positive finite number",
                   fixed = TRUE, label = nm)
    }
  }
})

test_that("the SCAD and MCP step error prints the scaled step", {
  msg_of <- function(expr) {
    tryCatch({
      expr
      ""
    }, error = function(e) conditionMessage(e))
  }
  # under the map d = (0.5, 0.25) the scaled steps are t d^2, and the largest
  # is 10.9 * 0.25 = 2.725 for SCAD and 12.4 * 0.25 = 3.1 for MCP
  D <- Matrix::Diagonal(x = c(0.5, 0.25))
  m <- msg_of(penalty_prox(scad_penalty(map = D), c(1, 1), 10.9,
                           list(lambda = 1, a = 3.7)))
  expect_match(m, "step < a - 1", fixed = TRUE)
  expect_match(m, "2.725", fixed = TRUE)
  m <- msg_of(penalty_prox(mcp_penalty(map = D), c(1, 1), 12.4,
                           list(lambda = 1, gamma = 3)))
  expect_match(m, "step < gamma", fixed = TRUE)
  expect_match(m, "\\b3\\.1\\b")
  # with curvature = 2 the scaled step is 6 * 2 * 0.5^2 = 3
  m <- msg_of(penalty_prox(scad_penalty(map = Matrix::Diagonal(x = 0.5),
                                        curvature = 2), 1, 6,
                           list(lambda = 1, a = 3.7)))
  expect_match(m, "step < a - 1", fixed = TRUE)
  expect_match(m, "scaled step is 3\\b")
})

test_that("check_curv() writes the singular for one coordinate", {
  expect_error(scad_penalty(n_coef = 1, curvature = c(1, 2)),
               "'curvature' must have 1 value, one per coordinate; it has 2.",
               fixed = TRUE)
  expect_error(scad_penalty(n_coef = 3, curvature = c(1, 2)),
               "'curvature' must have 3 values, one per coordinate; it has 2.",
               fixed = TRUE)
})
