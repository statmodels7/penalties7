#' @include prox.R
NULL

#' The Proximal Operator as a Piecewise Linear Table
#'
#' @description
#' Returns the scalar proximal operator of a separable penalty as an odd
#' piecewise linear map: three matrices of cuts, slopes and intercepts from
#' which a compiled loop can apply the operator without knowing which family it
#' came from. Returns `NULL` for a penalty whose operator has no such
#' description, and a caller that gets `NULL` uses [penalty_prox()] itself.
#'
#' @details
#' # What the table says
#'
#' The operator of a separable penalty acts one coordinate at a time and is
#' odd, so it is determined by what it does to \eqn{\lvert u\rvert}: on the
#' \eqn{k}-th interval, \eqn{\lvert u\rvert \le} `cut[j, k]`,
#'
#' \deqn{\mathrm{prox}(u) = \mathrm{sign}(u)\,(a_{jk}\lvert u\rvert + b_{jk}).}
#'
#' Every closed form the package carries has that shape. At
#' \eqn{\lambda = 1.5} for the first three, \eqn{\lambda = 1}, \eqn{a = 3.7},
#' \eqn{\gamma = 3} for the last two, and a step of `0.5`:
#'
#' | penalty | pieces | cuts |
#' |---|---|---|
#' | Gaussian prior (ridge) | 1 | `Inf` |
#' | Laplace prior (lasso) | 2 | `0.75`, `Inf` |
#' | elastic net | 2 | `0.45`, `Inf` |
#' | MCP | 3 | `0.5`, `3`, `Inf` |
#' | SCAD | 4 | `0.5`, `1.5`, `3.7`, `Inf` |
#'
#' The last cut is always `Inf`. Applied with [prox_apply()] and compared
#' against [penalty_prox()] at every cut, on both sides of every cut, and over
#' a grid across the whole range, all five agree to `8.9e-16` or better.
#'
#' # Why a table rather than the operator
#'
#' A coordinate descent applies the operator once per coordinate per sweep, at
#' a point that moves every time, so a compiled loop calling back into R for it
#' would spend its gain on the calls. Passing the numbers lets the loop stay
#' compiled while the penalty keeps the mathematics: the kernel evaluates any
#' map of this shape and names no family.
#'
#' # Why the step is a vector
#'
#' In a coordinate descent the step of coordinate \eqn{j} is \eqn{1/v_j} with
#' \eqn{v_j = \sum_i w_i x_{ij}^2}, which does not move while the working
#' weights are held. The whole table is therefore built once per weighted least
#' squares iteration and every sweep reads it. Note the asymmetry with
#' [penalty_prox()], which takes a **single** step.
#'
#' # Under a diagonal map
#'
#' Standardization is a diagonal \eqn{D}, under which a separable penalty stays
#' separable and the table survives. The builder is called at the step
#' \eqn{t_j d_j^2}, and the resulting cuts and intercepts are divided by
#' \eqn{\lvert d_j\rvert} while the slopes do not move, the slope multiplying a
#' point that was scaled and then divided back. The operator is odd, so only
#' the magnitude of \eqn{d} enters. For a lasso this leaves the cut at
#' \eqn{t_j d_j \lambda}: at \eqn{d = (0.5, 2, 3)}, \eqn{t = 0.4} and
#' \eqn{\lambda = 1.2} the cuts are `0.24`, `0.96`, `1.44`.
#'
#' The convexity condition of SCAD and MCP is tested on the scaled step, so it
#' tightens to \eqn{t < (a-1)/d_j^2} and \eqn{t < \gamma/d_j^2}: a standardized
#' penalty takes shorter steps.
#'
#' # What has no table
#'
#' `NULL` comes back for a quadratic or structured penalty, whose operator is
#' one linear solve over every coordinate at once; for a separable penalty under
#' a
#' map that is not diagonal, which is the generalized-lasso problem; for one
#' whose parent is read blockwise, whose coordinates do not separate; for one
#' whose operator is a root rather than a formula, such as the heavy-tailed
#' Student t prior; for one whose parent is not centered where the quadratic
#' pull is; and for SCAD or MCP at a step where the operator is set-valued,
#' that is `step >= a - 1` or `step >= gamma`.
#'
#' @param pen A [penalty()] object of any branch.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param step A numeric vector of step lengths, one per coefficient, recycled
#'   from a single value. Every entry must be positive, and **that is not
#'   checked**: a negative step produces a negative cut and a table nothing
#'   will reject. [penalty_prox()], which takes a single step, does check it.
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
#' # The table reproduces the operator exactly. The first point sits on its
#' # own cut and is set to zero; the second survives, shrunk by t * lambda.
#' u <- c(0.75, 4)
#' prox_apply(sp, u)
#' c(penalty_prox(lasso_penalty(n_coef = 1), u[1], 0.5, list(lambda = 1.5)),
#'   penalty_prox(lasso_penalty(n_coef = 1), u[2], 2, list(lambda = 1.5)))
#'
#' # SCAD needs four pieces, and has none where its operator is set-valued.
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
#' Builds the three matrices [penalty_prox_spec()] returns from a function that
#' gives one column per piece at a single step, recycling a step that is the
#' same for every coefficient.
#'
#' @details
#' Every branch's builder differs only in that per-step function, so this holds
#' the recycling and the transposition once. The number of pieces is read from
#' the first coefficient's answer and every other is required to match, which
#' is automatic: the piece count depends on the family and on the
#' hyperparameters, never on the step.
#'
#' @param step The step lengths, a numeric vector recycled to `n_coef`.
#' @param n_coef How many coefficients, a single whole number.
#' @param pieces A function of one step returning a matrix of three rows,
#'   `cut`, `slope` and `icept`, and one column per piece.
#'
#' @return A list of three matrices, `cut`, `slope` and `icept`, each `n_coef`
#'   by the piece count.
#'
#' @seealso [penalty_prox_spec()], [prox_apply()]
#'
#' @keywords internal
prox_table <- function(step, n_coef, pieces) {
  step <- rep_len(as.numeric(step), n_coef)
  rows <- lapply(step, pieces)
  k <- ncol(rows[[1L]])
  m <- function(i) matrix(vapply(rows, function(r) r[i, ], numeric(k)),
                          nrow = n_coef, ncol = k, byrow = TRUE)
  list(cut = m(1L), slope = m(2L), icept = m(3L))
}

#' Carry a Table Through a Diagonal Map
#'
#' @description
#' Builds the table of a separable penalty under a diagonal map from the
#' builder of its identity-map table, and returns `NULL` where there is a map
#' and it is not diagonal. With no map at all the builder is called directly.
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
#' The step handed to the builder is \eqn{t d^2}, which is where the convexity
#' condition of SCAD and MCP tightens: a standardized penalty takes shorter
#' steps and a table may come back `NULL` where the unmapped one exists.
#'
#' The map is recognized by its class through [map_diagonal()], so a base
#' matrix that happens to be diagonal is not treated as one.
#'
#' @param pen A [penalty()] object whose map is `NULL` or diagonal.
#' @param step A numeric vector of step lengths, recycled to `pen@n_coef`.
#' @param build A function of a penalty and a step returning the table, or
#'   `NULL`. Called on the penalty with its map removed.
#'
#' @return The list [penalty_prox_spec()] returns, or `NULL` when the map is
#'   not diagonal or the builder itself declines.
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
  lam <- theta$lambda
  a <- theta$a
  if (any(step >= a - 1)) return(NULL)
  prox_table(step, as.integer(pen@n_coef), function(t) {
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
  lam <- theta$lambda
  gam <- theta$gamma
  if (any(step >= gam)) return(NULL)
  prox_table(step, as.integer(pen@n_coef), function(t) {
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
#' Evaluates the map [penalty_prox_spec()] describes, in R. This is what a
#' compiled coordinate descent does to each coordinate, and what the tests
#' compare the table against.
#'
#' @details
#' For coordinate \eqn{j} it finds the first cut that \eqn{\lvert u_j\rvert}
#' does not exceed and returns
#' \eqn{\mathrm{sign}(u_j)(a_{jk}\lvert u_j\rvert + b_{jk})} for that piece.
#' The last cut is `Inf`, so a piece is always found; the `NA` guard covers a
#' table whose cuts are not increasing, which no shipped branch produces.
#'
#' Comparing this against [penalty_prox()] is how a table is verified. Over the
#' five families that have one, probed at every cut, on both sides of every cut
#' and across the whole range, the worst disagreement is `8.9e-16`.
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
