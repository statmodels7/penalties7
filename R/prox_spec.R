#' @include prox.R
NULL

#' The Proximal Operator as a Piecewise Linear Table
#'
#' @description
#' Describes the scalar proximal operator of a separable penalty as an odd
#' piecewise linear map, so that a compiled loop can apply it without knowing
#' which family it came from.
#'
#' @details
#' \strong{What the table says.} The operator of a separable penalty acts one
#' coordinate at a time and is odd, so it is determined by what it does to
#' \eqn{|u|}: on the \eqn{k}-th interval, \eqn{|u| \le} \code{cut[j, k]},
#' \deqn{\mathrm{prox}(u) = \mathrm{sign}(u)\,(a_{jk}|u| + b_{jk}).}
#' Every closed form the package carries has this shape: the soft threshold is
#' two pieces, the elastic net two, MCP three and SCAD four, and a Gaussian
#' prior is the single piece \eqn{u/(1 + t/\sigma^2)}.
#'
#' \strong{Why a table rather than the operator.} A coordinate descent applies
#' the operator once per coordinate per sweep, at a point that moves every
#' time, so a compiled loop calling back into \R for it would spend its gain on
#' the calls -- the property that decided how far the score-driven filter of
#' \pkg{modelterms7} could be ported. Passing the numbers instead lets the loop
#' stay compiled while the penalty keeps the mathematics: the kernel evaluates
#' any map of this shape and names no family.
#'
#' \strong{Why the step is a vector.} In a coordinate descent the step of
#' coordinate \eqn{j} is \eqn{1/v_j} with \eqn{v_j = \sum_i w_i x_{ij}^2},
#' which does not move while the working weights are held, so the whole table
#' is built once per weighted least squares iteration and the sweeps read it.
#'
#' \strong{A diagonal map.} Standardization is a diagonal \eqn{D}, under which
#' a separable penalty stays separable and the table survives: reading it at
#' \eqn{d_j v_j} with the step \eqn{t_j d_j^2} and dividing back gives the
#' cuts and the intercepts divided by \eqn{|d_j|} and the slopes unchanged.
#' The convexity condition of SCAD and MCP is tested on the scaled step, so it
#' becomes \eqn{t < (a-1)/d_j^2} and \eqn{t < \gamma/d_j^2}.
#'
#' \strong{What has no table.} A quadratic penalty under a general matrix is
#' not separable and returns \code{NULL}, as does a separable penalty under a
#' map that is not diagonal, one whose operator is a root rather than a
#' formula, and one whose parent is not centered where the quadratic pull is.
#' A caller that gets \code{NULL} uses \code{\link{penalty_prox}} itself.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param theta A named list of hyperparameter values.
#' @param step A numeric vector of step lengths, one per coefficient.
#' @param ... Passed to methods.
#'
#' @return A list with \code{cut}, \code{slope} and \code{icept}, each a matrix
#'   with one row per coefficient and one column per piece, the last cut being
#'   \code{Inf}; or \code{NULL} where the operator has no such description.
#'
#' @examples
#' pen <- lasso_penalty(n_coef = 2L)
#' penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))
#'
#' @seealso \code{\link{penalty_prox}}, \code{\link{has_prox}}
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
#' Builds the three matrices from one row per piece, recycling a step that is
#' the same for every coefficient.
#'
#' @param step The step lengths.
#' @param n_coef How many coefficients.
#' @param pieces A function of one step returning a matrix whose rows are
#'   \code{cut}, \code{slope}, \code{icept}.
#'
#' @return The list \code{\link{penalty_prox_spec}} returns.
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
#' builder of its identity-map table, and returns \code{NULL} where the map
#' is not diagonal.
#'
#' @details
#' A diagonal map only rescales each coordinate, and the identity
#' \deqn{\mathrm{prox}_{t\rho(d\,\cdot)}(v) =
#'   \mathrm{prox}_{t d^2 \rho}(d v)/d}
#' carries the table across. The identity-map table reads
#' \eqn{|w| \le \mathrm{cut}} to
#' \eqn{\mathrm{sign}(w)(\mathrm{slope}\,|w| + \mathrm{icept})} at
#' \eqn{w = dv}, so a cut on \eqn{|w|} is a cut on \eqn{|v|} divided by
#' \eqn{|d|} and the intercept divides by the same, while the slope, which
#' multiplies a point that was scaled and then divided back, does not move.
#' The operator is odd, so only the magnitude of \eqn{d} enters.
#'
#' The step handed to the builder is \eqn{t d^2}, which is where the
#' convexity condition of SCAD and MCP tightens: a standardized penalty
#' takes shorter steps.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param step A numeric vector of step lengths.
#' @param build A function of a penalty and a step returning the table, or
#'   \code{NULL}.
#'
#' @return The list \code{\link{penalty_prox_spec}} returns, or \code{NULL}.
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
  spec_diag(pen, step, function(pen, step) {
  fam <- .prox_family(pen)
  q <- as.integer(pen@n_coef)
  # the closed forms are written for a parent centered where the quadratic
  # pull is, and the table is odd, so an off-centre parent has none
  g0 <- penalty_gradient(pen, rep(0, q), theta)
  if (identical(fam, "gaussian1")) {
    if (max(abs(g0)) > 1e-8) return(NULL)
    s <- theta[[which(pen@params == "sigma")]]
    return(prox_table(step, q, function(t)
      rbind(Inf, 1 / (1 + t / s^2), 0)))
  }
  if (identical(fam, "laplace2") || identical(fam, "laplace")) {
    lam <- if (identical(fam, "laplace2")) theta[[which(pen@params == "lambda")]]
      else 1 / theta[[which(pen@params == "sigma")]]
    if (max(abs(g0)) > 1e-8 * max(1, lam)) return(NULL)
    return(prox_table(step, q, function(t)
      rbind(c(t * lam, Inf), c(0, 1), c(0, -t * lam))))
  }
  if (identical(fam, "enet")) {
    lam <- theta[[which(pen@params == "lambda")]]
    al <- theta[[which(pen@params == "alpha")]]
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
#' Evaluates the map \code{\link{penalty_prox_spec}} describes, in \R, which is
#' what a compiled loop does and what the tests compare the table against.
#'
#' @param spec A table, as \code{\link{penalty_prox_spec}} returns it.
#' @param u The points, one per coefficient.
#'
#' @return A numeric vector as long as \code{u}.
#'
#' @examples
#' pen <- lasso_penalty(n_coef = 2L)
#' sp <- penalty_prox_spec(pen, list(lambda = 1.5), step = c(0.5, 2))
#' prox_apply(sp, c(2, 2))
#'
#' @seealso \code{\link{penalty_prox_spec}}
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
