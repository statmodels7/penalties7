#' @include generics.R quadratic_penalty.R distrib_penalty.R scad_mcp.R
#' @include structured_penalty.R
NULL

#' @title Proximal Operator of a Penalty
#'
#' @description
#' The point that minimizes the penalty plus a quadratic pull towards
#' \eqn{v},
#' \deqn{\mathrm{prox}_{t\rho}(v) = \arg\min_{\beta}
#'   \left\{ \tfrac{1}{2t}\lVert \beta - v \rVert^{2} +
#'   \rho(\beta;\theta) \right\},}
#' which is what a proximal-gradient method applies after each gradient
#' step, and the only operation that lets a non-differentiable penalty be
#' minimized without differencing it.
#'
#' @details
#' Where the operator has a closed form, that form is used. A quadratic
#' penalty gives one linear solve, \eqn{(I + tS)^{-1}v} with \eqn{S} the
#' penalty's Hessian, and works for any linear map. The separable branches
#' need the identity map, because with a general \eqn{D} the problem does
#' not split by coordinate and has no elementary solution: that case is
#' rejected rather than approximated.
#'
#' The named separable instances are exact: the Gaussian gives
#' \eqn{v/(1+t/\sigma^{2})} and the Laplace the soft threshold
#' \eqn{\mathrm{sign}(v)(\lvert v \rvert - t\lambda)_{+}}. SCAD and MCP
#' have the closed piecewise operators of their defining papers, each
#' valid only while the quadratic pull dominates the concave region of the
#' penalty -- \eqn{t < a-1} for SCAD and \eqn{t < \gamma} for MCP -- so a
#' step beyond that is rejected, the operator being set-valued there.
#'
#' Any other separable penalty built from a distribution with a
#' differentiable log-density is solved coordinatewise from the stationary
#' condition \eqn{(\beta - v)/t = \ell^{(y)}(\beta)}. The right-hand side
#' is the response derivative
#' \code{\link[distributions7]{distrib_grad_y}}, closed form for every
#' continuous family, and log-concavity of the density makes the left side
#' minus the right side strictly increasing, so the root is unique and a
#' bracketed search finds it.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param v A numeric vector, the point to be pulled towards.
#' @param step The positive step length \eqn{t}.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of the same length as \code{v}.
#'
#' @seealso \code{\link{has_prox}}, \code{\link{penalty_value}}
#'
#' @examples
#' # the lasso: a soft threshold
#' penalty_prox(lasso_penalty(n_coef = 3), c(2, 0.3, -1.4), 1, list(lambda = 1))
#'
#' # the ridge: a shrinkage
#' penalty_prox(ridge_penalty(n_coef = 2), c(2, -1), 1, list(sigma = 1))
#'
#' @export
penalty_prox <- S7::new_generic("penalty_prox", "pen",
  function(pen, v, step, theta, ...) {
    theta <- align_ptheta(pen, theta)
    v <- as.numeric(v)
    if (length(step) != 1L || is.na(step) || step <= 0) {
      stop("'step' must be a single positive number.", call. = FALSE)
    }
    if (length(v) != pen@n_coef) {
      stop(sprintf("'v' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })

S7::method(penalty_prox, penalty) <- function(pen, v, step, theta, ...) {
  stop(sprintf(paste0(
    "'%s' has no proximal operator.\n",
    "  A penalty supplies one when it is quadratic, or separable under the\n",
    "  identity map. Use has_prox() to ask before calling."),
    pen@penalty_name), call. = FALSE)
}

#' @title Does a Penalty Supply a Proximal Operator?
#'
#' @description
#' \code{TRUE} when \code{\link{penalty_prox}} can be evaluated for this
#' penalty at a suitable step, so that a caller may choose a proximal
#' method over a smooth one without provoking an error.
#'
#' @details
#' The operator in question is
#'
#' \deqn{\operatorname{prox}_{t\rho}(v)
#'   = \arg\min_{\beta} \Bigl\{ \tfrac{1}{2}\lVert \beta - v \rVert^{2}
#'     + t\, \rho(\beta; \theta) \Bigr\},}
#'
#' which a proximal gradient method evaluates once per iteration and which
#' therefore has to be available in closed form, or as a solve, for the
#' method to be worth using. A penalty carries one when it is quadratic or
#' structured (one linear solve at any map), when it is separable with a
#' parent whose operator is closed or whose stationarity condition has a
#' coordinatewise root, or when it is SCAD or MCP over their convex
#' regions.
#'
#' @param pen A \code{\link{penalty}} object.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{penalty_prox}}
#'
#' @examples
#' c(lasso = has_prox(lasso_penalty()), scad = has_prox(scad_penalty()))
#'
#' @export
has_prox <- function(pen) {
  if (!S7::S7_inherits(pen, penalty)) {
    stop("'pen' must be a penalty object.", call. = FALSE)
  }
  m <- S7::method(penalty_prox, S7::S7_class(pen))
  owner <- attr(attr(m, "signature")[[1L]], "name")
  if (identical(owner, "penalty")) return(FALSE)
  if (is_quadratic(pen) || S7::S7_inherits(pen, StructuredPenalty)) {
    return(TRUE)
  }
  is.null(pen@map) || !is.null(map_diagonal(pen))
}

#' The Diagonal of a Map That Has One
#'
#' @description
#' The entries of \eqn{D} where the map is diagonal, and \code{NULL} where
#' there is no map or the map is not diagonal.
#'
#' @details
#' A diagonal map RESCALES each coordinate on its own, and a separable
#' penalty under one is still separable. That is what standardization comes
#' to: penalizing a column divided by its own spread is penalizing
#' \eqn{\rho(s_j\beta_j)}, so the scaling never has to touch the design and
#' the sparsity of a block survives it. A general \eqn{D} mixes coordinates
#' and turns the problem into the generalized-lasso one, which is a different
#' algorithm rather than a different formula.
#'
#' The map is recognized by its CLASS and not by inspecting its entries: a
#' \pkg{Matrix} diagonal object says what it is and costs \eqn{q} numbers,
#' where testing a dense matrix for diagonality would cost \eqn{q^2} and
#' defeat the point.
#'
#' @param pen A \code{\link{penalty}} object.
#'
#' @return A numeric vector, or \code{NULL}.
#'
#' @keywords internal
map_diagonal <- function(pen) {
  m <- pen@map
  if (is.null(m)) return(NULL)
  if (!isS4(m) || !methods::is(m, "diagonalMatrix")) return(NULL)
  d <- if (identical(m@diag, "U")) rep(1, nrow(m)) else as.numeric(m@x)
  if (anyNA(d) || any(d == 0)) return(NULL)
  d
}

# The separable branches split by coordinate under the identity map and under
# a DIAGONAL one, which only rescales each coordinate: with u = d b,
#
#     argmin_b (b - v)^2/(2t) + rho(d b)
#       = argmin_u (u - d v)^2/(2 t d^2) + rho(u),   b = u/d
#
# so the operator is the identity-map one read at the scaled point with the
# step scaled by d^2.
# the same penalty with its map removed, which is what the scaled point is
# handed to: the diagonal has already been applied to the point and the step
.undiag <- function(pen) {
  pen@map <- NULL
  pen
}

.prox_scaling <- function(pen) {
  if (is.null(pen@map)) return(NULL)
  d <- map_diagonal(pen)
  if (is.null(d)) {
    stop(sprintf(paste0(
      "'%s' carries a linear map that is not diagonal, and its proximal\n",
      "  operator does not split by coordinate. Only the identity map and a\n",
      "  diagonal one have a closed form here."),
      pen@penalty_name), call. = FALSE)
  }
  d
}

#' @title Proximal Operator of a Quadratic Penalty
#' @name penalty_prox.QuadraticPenalty
#' @description
#' One linear solve, \eqn{(I + tS)^{-1}v} with \eqn{S = \lambda D'PD} the
#' penalty's Hessian, which holds for any map because the objective stays
#' quadratic.
#' @param pen A \code{QuadraticPenalty} object.
#' @param v A numeric vector.
#' @param step The step length.
#' @param theta A named list of hyperparameter values.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(penalty_prox, QuadraticPenalty) <- function(pen, v, step, theta, ...) {
  S <- penalty_hessian(pen, v, theta)
  as.numeric(solve(diag(length(v)) + step * S, v))
}

#' @rdname penalty_prox.QuadraticPenalty
#' @name penalty_prox.StructuredPenalty
#' @keywords internal
S7::method(penalty_prox, StructuredPenalty) <- function(pen, v, step, theta, ...) {
  S <- penalty_hessian(pen, v, theta)
  as.numeric(solve(diag(length(v)) + step * S, v))
}

# The parent of a separable penalty, with any fixed() wrapper unwrapped, so
# that the closed forms below can be recognized by the family underneath.
.prox_family <- function(pen) {
  d <- pen@parent
  nm <- d@distrib_name
  sub("^fixed ", "", sub(" \\[.*$", "", nm))
}

#' @title Proximal Operator of a Separable Penalty
#' @name penalty_prox.DistribPenalty
#' @description
#' Closed form for the Gaussian (a shrinkage), the Laplace (a soft
#' threshold) and the elastic net (the one followed by the other); for any
#' other parent with a differentiable log-density, the
#' coordinatewise root of \eqn{(\beta - v)/t = \ell^{(y)}(\beta)}, which is
#' unique because the density is log-concave.
#' @param pen A \code{DistribPenalty} object.
#' @param v A numeric vector.
#' @param step The step length.
#' @param theta A named list of hyperparameter values.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(penalty_prox, DistribPenalty) <- function(pen, v, step, theta, ...) {
  d <- .prox_scaling(pen)
  if (!is.null(d)) {
    return(S7::method(penalty_prox, DistribPenalty)(
      .undiag(pen), d * v, step * d^2, theta, ...) / d)
  }
  fam <- .prox_family(pen)

  # The closed forms below hold for a parent centered where the quadratic
  # pull is, and a family name does not say where that is. The gradient at
  # the origin does: it is the location term of the stationary condition, so
  # the Gaussian carries it through, and the Laplace, whose location cannot
  # be recovered from one evaluation, is required to sit at zero.
  g0 <- penalty_gradient(pen, rep(0, length(v)), theta)

  if (identical(fam, "gaussian1")) {
    s <- theta[[which(pen@params == "sigma")]]
    return((v - step * g0) / (1 + step / s^2))
  }
  if (identical(fam, "laplace2") || identical(fam, "laplace")) {
    lam <- if (identical(fam, "laplace2")) {
      theta[[which(pen@params == "lambda")]]
    } else {
      1 / theta[[which(pen@params == "sigma")]]
    }
    if (max(abs(g0)) > 1e-8 * max(1, lam)) {
      stop(paste("the Laplace proximal operator is written for a parent",
                 "centered at zero; this one is not."), call. = FALSE)
    }
    return(sign(v) * pmax(abs(v) - step * lam, 0))
  }
  if (identical(fam, "enet")) {
    lam <- theta[[which(pen@params == "lambda")]]
    al <- theta[[which(pen@params == "alpha")]]
    if (max(abs(g0)) > 1e-8 * max(1, lam)) {
      stop(paste("the elastic-net proximal operator is written for a parent",
                 "centered at zero; this one is not."), call. = FALSE)
    }
    # the stationary condition (b - v)/t + a sgn(b) + c b = 0 separates
    # into the soft threshold of the Laplace part followed by the
    # shrinkage of the Gaussian one
    sft <- sign(v) * pmax(abs(v) - step * lam * al, 0)
    return(sft / (1 + step * lam * (1 - al)))
  }

  if (length(pen@kinks)) {
    stop(sprintf(paste0(
      "'%s' declares a kink and is not one of the closed-form families, so\n",
      "  its proximal operator cannot be found by a smooth root."),
      pen@penalty_name), call. = FALSE)
  }

  # h(b) = (b - v)/t - l_y(b) is strictly increasing for a log-concave
  # density, so one bracketed root per coordinate.
  ly <- function(b) {
    as.numeric(distributions7::distrib_grad_y(pen@parent, b, theta)) + 0 * b
  }
  vapply(v, function(vi) {
    h <- function(b) (b - vi) / step - ly(b)
    lo <- min(vi, 0) - 1
    hi <- max(vi, 0) + 1
    for (k in 1:60) {
      if (h(lo) <= 0 && h(hi) >= 0) break
      if (h(lo) > 0) lo <- lo - (hi - lo)
      if (h(hi) < 0) hi <- hi + (hi - lo)
    }
    if (h(lo) > 0 || h(hi) < 0) {
      stop("the proximal root could not be bracketed; is the density log-concave?",
           call. = FALSE)
    }
    stats::uniroot(h, c(lo, hi), tol = .Machine$double.eps^0.75)$root
  }, numeric(1))
}

#' @title Proximal Operator of SCAD
#' @name penalty_prox.ScadPenalty
#' @description
#' The closed piecewise operator: a soft threshold near zero, a rescaled
#' threshold on the tapering region, and the identity beyond \eqn{a\lambda}.
#' The middle region is convex only while \eqn{t < a - 1}, and a longer
#' step is rejected.
#' @param pen A \code{ScadPenalty} object.
#' @param v A numeric vector.
#' @param step The step length.
#' @param theta A named list containing \code{lambda} and \code{a}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(penalty_prox, ScadPenalty) <- function(pen, v, step, theta, ...) {
  d <- .prox_scaling(pen)
  if (!is.null(d)) {
    return(S7::method(penalty_prox, ScadPenalty)(
      .undiag(pen), d * v, step * d^2, theta, ...) / d)
  }
  lam <- theta$lambda
  a <- theta$a
  if (any(step >= a - 1)) {
    stop(sprintf(paste0(
      "the SCAD proximal operator needs step < a - 1 (%g here): beyond that the\n",
      "  subproblem is not convex and the operator is set-valued."), a - 1),
      call. = FALSE)
  }
  s <- sign(v)
  u <- abs(v)
  ifelse(u <= (1 + step) * lam,
         s * pmax(u - step * lam, 0),
         ifelse(u <= a * lam,
                s * (u - step * a * lam / (a - 1)) / (1 - step / (a - 1)),
                v))
}

#' @title Proximal Operator of MCP
#' @name penalty_prox.McpPenalty
#' @description
#' The closed piecewise operator: a rescaled soft threshold below
#' \eqn{\gamma\lambda} and the identity beyond it. Convex only while
#' \eqn{t < \gamma}, and a longer step is rejected.
#' @param pen An \code{McpPenalty} object.
#' @param v A numeric vector.
#' @param step The step length.
#' @param theta A named list containing \code{lambda} and \code{gamma}.
#' @param ... Unused.
#' @return A numeric vector.
#' @keywords internal
S7::method(penalty_prox, McpPenalty) <- function(pen, v, step, theta, ...) {
  d <- .prox_scaling(pen)
  if (!is.null(d)) {
    return(S7::method(penalty_prox, McpPenalty)(
      .undiag(pen), d * v, step * d^2, theta, ...) / d)
  }
  lam <- theta$lambda
  gam <- theta$gamma
  if (any(step >= gam)) {
    stop(sprintf(paste0(
      "the MCP proximal operator needs step < gamma (%g here): beyond that the\n",
      "  subproblem is not convex and the operator is set-valued."), gam),
      call. = FALSE)
  }
  s <- sign(v)
  u <- abs(v)
  ifelse(u <= gam * lam,
         s * pmax(u - step * lam, 0) / (1 - step / gam),
         v)
}
