#' @include prox.R
NULL

#' The Proximal Operator as a Piecewise Linear Table
#'
#' @description
#' Returns the scalar proximal operator of a separable penalty as an odd
#' piecewise linear map: three matrices of cuts, slopes and intercepts from
#' which a compiled loop can apply the operator without reference to the
#' family. Returns `NULL` for a penalty whose operator has no such
#' description; where such a penalty has an operator, it is applied with
#' [penalty_prox()].
#'
#' @details
#' # The piecewise linear map
#'
#' The operator of a separable penalty acts one coordinate at a time and is
#' odd, so it is determined by its action on \eqn{\lvert u\rvert}: on the
#' \eqn{k}-th interval, the first for which \eqn{\lvert u\rvert \le}
#' `cut[j, k]`,
#'
#' \deqn{\mathrm{prox}(u) = \mathrm{sign}(u)\,(a_{jk}\lvert u\rvert + b_{jk}),}
#'
#' with \eqn{a_{jk}} = `slope[j, k]` and \eqn{b_{jk}} = `icept[j, k]`. Every
#' closed-form operator of a separable penalty in the package has this shape.
#' The table below is for a step of `0.5`, with \eqn{\sigma = 0.8} for the
#' Gaussian prior (whose single cut does not depend on \eqn{\sigma}),
#' \eqn{\lambda = 1.5} for the Laplace prior, \eqn{\lambda = 1.5} and
#' \eqn{\alpha = 0.6} for the elastic net, \eqn{\lambda = 1} and
#' \eqn{\gamma = 3} for MCP, and \eqn{\lambda = 1} and \eqn{a = 3.7} for SCAD:
#'
#' | penalty | pieces | cuts |
#' |---|---|---|
#' | Gaussian prior (`gaussian1` parent) | 1 | `Inf` |
#' | Laplace prior ([lasso_penalty()]) | 2 | `0.75`, `Inf` |
#' | elastic net ([elasticnet_penalty()]) | 2 | `0.45`, `Inf` |
#' | MCP | 3 | `0.5`, `3`, `Inf` |
#' | SCAD | 4 | `0.5`, `1.5`, `3.7`, `Inf` |
#'
#' The last cut is always `Inf`. [prox_apply()] evaluates a table, and the
#' tests of the package compare it with [penalty_prox()] at every cut, on both
#' sides of every cut and over a grid of points.
#'
#' # Accuracy near the step limit of SCAD and MCP
#'
#' The table of SCAD and MCP is ill-conditioned as the step approaches its
#' limit, \eqn{a - 1} for SCAD and \eqn{\gamma} for MCP. On the tapering piece
#' the slope is \eqn{1/(1 - t/(a-1))} for SCAD and \eqn{1/(1 - t/\gamma)} for
#' MCP, which grows without bound as the step approaches the limit, and the
#' intercept grows with it and has the opposite sign, so
#' \eqn{a_{jk}\lvert u\rvert + b_{jk}} is the difference of two large numbers
#' that nearly cancel. The table is therefore less accurate there than
#' [penalty_prox()], which evaluates the same piece in a form without this
#' cancellation. With a curvature or a diagonal map the limit applies to the
#' scaled step \eqn{t_j c_j d_j^2} described below.
#'
#' # Why a table rather than the operator
#'
#' A coordinate descent applies the operator once per coordinate per sweep, at
#' a point that moves every time, so a compiled loop that called back into R
#' for it would lose its gain to the calls. Passing the numbers keeps the loop
#' compiled. The table holds the mathematics of the penalty, and the compiled
#' kernel evaluates any map of this shape without reference to the family.
#'
#' # Why the step is a vector
#'
#' In a coordinate descent the step of coordinate \eqn{j} is \eqn{1/v_j} with
#' \eqn{v_j = \sum_i w_i x_{ij}^2}, which does not change while the working
#' weights are held fixed. The whole table is therefore built once per
#' weighted least squares iteration and read by every sweep. [penalty_prox()],
#' by contrast, takes a **single** step.
#'
#' # Under a diagonal map
#'
#' Standardization is a diagonal \eqn{D}, under which a separable penalty stays
#' separable and the table is still available. The builder is called at the
#' step \eqn{t_j d_j^2}, and the cuts and intercepts that it returns are
#' divided by \eqn{\lvert d_j\rvert}, while the slopes are unchanged because
#' each slope multiplies a point that was scaled and then divided back. The
#' operator is odd, so only the magnitude of \eqn{d} enters. For a lasso the
#' cut is \eqn{t_j \lvert d_j\rvert \lambda} (with \eqn{d = (0.5, 2, 3)},
#' \eqn{t = 0.4} and \eqn{\lambda = 1.2} the cuts are `0.24`, `0.96` and
#' `1.44`).
#'
#' The step limit of SCAD and MCP is tested on the scaled step. With
#' \eqn{c_j} the curvature of the coordinate (the `curv` property, one when it
#' is empty), the table exists only while \eqn{t_j c_j d_j^2 < a - 1} for SCAD
#' and \eqn{t_j c_j d_j^2 < \gamma} for MCP in every coordinate, so a
#' standardized penalty takes shorter steps.
#'
#' A diagonal map is recognized by its class through [map_diagonal()].
#' [distrib_penalty()] (and therefore [lasso_penalty()] and
#' [elasticnet_penalty()]), [scad_penalty()] and [mcp_penalty()] store a
#' square map whose off-diagonal entries are all zero as a
#' [Matrix::Diagonal()] object, so a diagonal base matrix gives a table. A
#' diagonal map with a zero entry gives `NULL`.
#'
#' # Penalties without a table
#'
#' `NULL` is returned for a quadratic or structured penalty, whose operator is
#' one linear solve over every coordinate at once, and for an additive
#' penalty. Among separable penalties, `NULL` is returned under a map that is
#' not recognized as diagonal, which is the generalized-lasso problem; for a
#' parent read blockwise, whose coordinates do not separate; for a parent
#' other than the Gaussian, Laplace and elastic-net families, whose operator
#' has no closed form (the Student t prior, for example, has an operator
#' computed by root finding); for a parent whose gradient at zero is not zero,
#' that is, a parent not centered at zero; and for SCAD or MCP when the scaled
#' step \eqn{t_j c_j d_j^2} of some coordinate reaches \eqn{a - 1} or
#' \eqn{\gamma}, where the operator is set-valued.
#'
#' @param pen A [penalty()] object of any branch.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param step A numeric vector of step lengths, either a single value used for
#'   every coefficient or one value per coefficient. Every entry must be finite
#'   and positive, and any other length or value is rejected with an error
#'   naming `step`.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A list of three matrices, `cut`, `slope` and `icept`, each with
#'   `pen@n_coef` rows and one column per piece, the last column of `cut` being
#'   `Inf`. `NULL` where the operator has no piecewise linear description.
#'
#' @examples
#' # A lasso at two different steps: the threshold is t * lambda.
#' pen <- lasso_penalty(n_coef = 2L)
#' sp <- penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))
#' sp$cut
#'
#' # The table reproduces the operator. The first point lies on its cut and
#' # is set to zero; the second is shrunk by t * lambda.
#' u <- c(0.75, 4)
#' prox_apply(sp, u)
#' c(penalty_prox(lasso_penalty(n_coef = 1), u[1], 0.5, list(lambda = 1.5)),
#'   penalty_prox(lasso_penalty(n_coef = 1), u[2], 2, list(lambda = 1.5)))
#'
#' # SCAD has four pieces, and no table at a step where its operator is
#' # set-valued.
#' penalty_prox_spec(scad_penalty(), list(lambda = 1, a = 3.7),
#'                   step = 0.5)$cut
#' penalty_prox_spec(scad_penalty(), list(lambda = 1, a = 3.7), step = 2.7)
#'
#' # A quadratic penalty has a solve, not a coordinatewise map.
#' penalty_prox_spec(quadratic_penalty(diag(2)), list(lambda = 1), step = 0.5)
#'
#' @seealso [prox_apply()] to evaluate the table, [penalty_prox()] for the
#'   operator itself, [has_prox()] for whether there is one at all,
#'   [spec_diag()] for the diagonal-map transport.
#'
#' @export
penalty_prox_spec <- S7::new_generic("penalty_prox_spec", "pen",
  function(pen, theta, step, ...) {
    theta <- align_ptheta(pen, theta)
    n <- as.integer(pen@n_coef)
    if (!is.numeric(step) || !(length(step) %in% c(1L, n)) ||
        anyNA(step) || any(!is.finite(step)) || any(step <= 0)) {
      stop(sprintf(paste0(
        "'step' must be finite and positive, a single value or %d values,\n",
        "  one per coefficient."), n), call. = FALSE)
    }
    S7::S7_dispatch()
  })

#' @name penalty_prox_spec.penalty
#' @rdname penalty_prox_spec
#' @keywords internal
S7::method(penalty_prox_spec, penalty) <- function(pen, theta, step, ...) {
  NULL
}


#' Assemble a Piecewise Linear Table
#'
#' @description
#' Builds the three matrices that [penalty_prox_spec()] returns from a
#' function that gives one column per piece for a single coefficient.
#' `prox_table()` takes a function of the step alone, for a family whose
#' pieces are the same in every coordinate; `prox_table_j()` takes a function
#' of the step and the coordinate index, for a family whose rate can differ
#' by coordinate, as the rate \eqn{\lambda/c_j} of SCAD and MCP does.
#'
#' @details
#' The step vector is recycled to `n_coef`, and the function is called once
#' per coefficient with that coefficient's own step. The builders of the
#' branches differ only in that function, so the recycling and the
#' transposition are written here once. The number of pieces is read from the
#' first coefficient's columns, and the others must have the same number,
#' which holds for every shipped builder because the piece count depends on
#' the family and not on the step or the coordinate.
#'
#' @param step The step lengths, a numeric vector recycled to `n_coef`.
#' @param n_coef How many coefficients, a single whole number.
#' @param pieces A function returning a matrix of three rows, `cut`, `slope`
#'   and `icept`, and one column per piece. For `prox_table()` it takes one
#'   argument, the step; for `prox_table_j()` it takes two, the step and the
#'   coordinate index \eqn{j}.
#'
#' @return A list of three matrices, `cut`, `slope` and `icept`, each `n_coef`
#'   by the piece count.
#'
#' @seealso [penalty_prox_spec()], [prox_apply()]
#'
#' @keywords internal
prox_table <- function(step, n_coef, pieces) {
  prox_table_j(step, n_coef, function(t, j) pieces(t))
}

#' @rdname prox_table
#' @keywords internal
prox_table_j <- function(step, n_coef, pieces) {
  step <- rep_len(as.numeric(step), n_coef)
  rows <- lapply(seq_len(n_coef), function(j) pieces(step[[j]], j))
  k <- ncol(rows[[1L]])
  m <- function(i) matrix(vapply(rows, function(r) r[i, ], numeric(k)),
                          nrow = n_coef, ncol = k, byrow = TRUE)
  list(cut = m(1L), slope = m(2L), icept = m(3L))
}

#' Carry a Table Through a Diagonal Map
#'
#' @description
#' Builds the table of a separable penalty under a diagonal map from the
#' builder of its identity-map table, and returns `NULL` where the penalty has
#' a map that is not recognized as diagonal. With no map the builder is
#' called directly.
#'
#' @details
#' A diagonal map only rescales each coordinate, and the identity
#'
#' \deqn{\mathrm{prox}_{t\rho(d\,\cdot)}(v) =
#'   \mathrm{prox}_{t d^2 \rho}(d v)/d}
#'
#' carries the table across. The identity-map table reads
#' \eqn{\lvert w\rvert \le \mathrm{cut}} to
#' \eqn{\mathrm{sign}(w)(\mathrm{slope}\,\lvert w\rvert + \mathrm{icept})} at
#' \eqn{w = dv}, so a cut on \eqn{\lvert w\rvert} is a cut on
#' \eqn{\lvert v\rvert} divided by \eqn{\lvert d\rvert}, the intercept divides
#' by the same, and the slope does not move, multiplying a point that was
#' scaled and then divided back. The operator is odd, so only the magnitude of
#' \eqn{d} enters.
#'
#' The step passed to the builder is \eqn{t d^2}, and the step limit of SCAD
#' and MCP is tested on it, so a standardized penalty takes shorter steps and
#' the table can be `NULL` at a step where the table without the map exists.
#'
#' The map is recognized by its class through [map_diagonal()], which returns
#' `NULL` for a base matrix and for a diagonal with a zero entry. The
#' constructors [distrib_penalty()], [scad_penalty()] and [mcp_penalty()] store
#' a diagonal base matrix as a [Matrix::Diagonal()] object, so its diagonal
#' map is recognized.
#'
#' @param pen A [penalty()] object whose map is `NULL` or diagonal.
#' @param step A numeric vector of step lengths, of length one or
#'   `pen@n_coef`, recycled to `pen@n_coef`.
#' @param build A function of a penalty and a step returning the table, or
#'   `NULL`. It is called on the penalty with its map removed.
#'
#' @return The list that [penalty_prox_spec()] returns, or `NULL` when the
#'   map is not recognized as diagonal or the builder returns `NULL`.
#'
#' @seealso [penalty_prox_spec()], [map_diagonal()], [prox_table()]
#'
#' @keywords internal
spec_diag <- function(pen, step, build) {
  if (is.null(pen@map)) return(build(pen, step))
  d <- map_diagonal(pen)
  if (is.null(d)) return(NULL)
  q <- as.integer(pen@n_coef)
  # the operator is odd, so the sign of d cancels between the point and the
  # division back and only the magnitude enters
  d <- abs(rep_len(d, q))
  sp <- build(.undiag(pen), rep_len(as.numeric(step), q) * d^2)
  if (is.null(sp)) return(NULL)
  sp$cut <- sp$cut / d
  sp$icept <- sp$icept / d
  sp
}


#' @name penalty_prox_spec.DistribPenalty
#' @rdname penalty_prox_spec
#' @keywords internal
S7::method(penalty_prox_spec, DistribPenalty) <- function(pen, theta, step,
                                                          ...) {
  # a blockwise parent has no operator, so it has no table either
  if (pen@block > 1L) return(NULL)
  spec_diag(pen, step, function(pen, step) {
  fam <- .prox_family(pen)
  q <- as.integer(pen@n_coef)
  # the closed forms are written for a parent centered where the quadratic
  # pull is, and the table is odd, so an off-center parent has none
  g0 <- penalty_gradient(pen, rep(0, q), theta)
  if (identical(fam, "gaussian1")) {
    if (max(abs(g0)) > 1e-8) return(NULL)
    s <- .prox_param(pen, theta, "sigma")
    return(prox_table(step, q, function(t)
      rbind(Inf, 1 / (1 + t / s^2), 0)))
  }
  if (identical(fam, "laplace2") || identical(fam, "laplace")) {
    lam <- if (identical(fam, "laplace2")) .prox_param(pen, theta, "lambda")
      else 1 / .prox_param(pen, theta, "sigma")
    if (max(abs(g0)) > 1e-8 * max(1, lam)) return(NULL)
    return(prox_table(step, q, function(t)
      rbind(c(t * lam, Inf), c(0, 1), c(0, -t * lam))))
  }
  if (identical(fam, "enet")) {
    lam <- .prox_param(pen, theta, "lambda")
    al <- .prox_param(pen, theta, "alpha")
    if (max(abs(g0)) > 1e-8 * max(1, lam)) return(NULL)
    a <- lam * al
    cc <- lam * (1 - al)
    return(prox_table(step, q, function(t) {
      d <- 1 + t * cc
      rbind(c(t * a, Inf), c(0, 1 / d), c(0, -t * a / d))
    }))
  }
  NULL
  })
}


#' @name penalty_prox_spec.ScadPenalty
#' @rdname penalty_prox_spec
#' @keywords internal
S7::method(penalty_prox_spec, ScadPenalty) <- function(pen, theta, step, ...) {
  spec_diag(pen, step, function(pen, step) {
  n <- as.integer(pen@n_coef)
  cv <- curv_of(pen, n)
  L <- theta$lambda / cv
  a <- theta$a
  step <- rep_len(as.numeric(step), n) * cv
  if (any(step >= a - 1)) return(NULL)
  prox_table_j(step, n, function(t, j) {
    lam <- L[[j]]
    d <- 1 - t / (a - 1)
    rbind(c(t * lam, (1 + t) * lam, a * lam, Inf),
          c(0, 1, 1 / d, 1),
          c(0, -t * lam, -t * a * lam / ((a - 1) * d), 0))
  })
  })
}


#' @name penalty_prox_spec.McpPenalty
#' @rdname penalty_prox_spec
#' @keywords internal
S7::method(penalty_prox_spec, McpPenalty) <- function(pen, theta, step, ...) {
  spec_diag(pen, step, function(pen, step) {
  n <- as.integer(pen@n_coef)
  cv <- curv_of(pen, n)
  L <- theta$lambda / cv
  gam <- theta$gamma
  step <- rep_len(as.numeric(step), n) * cv
  if (any(step >= gam)) return(NULL)
  prox_table_j(step, n, function(t, j) {
    lam <- L[[j]]
    d <- 1 - t / gam
    rbind(c(t * lam, gam * lam, Inf),
          c(0, 1 / d, 1),
          c(0, -t * lam / d, 0))
  })
  })
}


#' Apply a Piecewise Linear Table
#'
#' @description
#' Evaluates in R the map that [penalty_prox_spec()] describes. A compiled
#' coordinate descent applies the same map to each coordinate, and the tests
#' compare the result with [penalty_prox()].
#'
#' @details
#' For coordinate \eqn{j} the function finds the first cut that
#' \eqn{\lvert u_j\rvert} does not exceed and returns
#' \eqn{\mathrm{sign}(u_j)(a_{jk}\lvert u_j\rvert + b_{jk})} for that piece.
#' The last cut of every shipped table is `Inf`, so a piece is always found.
#' If \eqn{\lvert u_j\rvert} exceeds every cut of its row, which can happen
#' only in a table whose last cut is finite, the last piece is used.
#'
#' A table is verified by comparing this function with [penalty_prox()] at
#' every cut, on both sides of every cut and over a grid of points. For SCAD
#' and MCP at a step close to its limit the table is less accurate than
#' [penalty_prox()], as the page of [penalty_prox_spec()] explains.
#'
#' @param spec A table, as [penalty_prox_spec()] returns it. Its matrices must
#'   have at least as many rows as `u` is long.
#' @param u The points, one per coefficient, a numeric vector.
#'
#' @return A numeric vector as long as `u`.
#'
#' @examples
#' pen <- lasso_penalty(n_coef = 2L)
#' sp <- penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))
#'
#' # The first coordinate is thresholded at 0.75 and the second at 3.
#' prox_apply(sp, c(2, 2))
#' prox_apply(sp, c(0.5, 2.9))
#'
#' # Which is exactly what the operator gives, coordinate by coordinate.
#' one <- lasso_penalty(n_coef = 1)
#' c(penalty_prox(one, 2, 0.5, list(lambda = 1.5)),
#'   penalty_prox(one, 2, 2, list(lambda = 1.5)))
#'
#' # The map is odd, so negating the point negates the answer.
#' prox_apply(sp, c(-2, -2))
#'
#' @seealso [penalty_prox_spec()] for the table, [penalty_prox()] for the
#'   operator it stands in for.
#'
#' @export
prox_apply <- function(spec, u) {
  vapply(seq_along(u), function(j) {
    a <- abs(u[[j]])
    k <- which(a <= spec$cut[j, ])[1L]
    if (is.na(k)) k <- ncol(spec$cut)
    sign(u[[j]]) * (spec$slope[j, k] * a + spec$icept[j, k])
  }, numeric(1))
}
