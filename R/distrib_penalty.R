#' @include generics.R
NULL

#' @title S7 Class for the Separable Penalty
#'
#' @description
#' The class \code{\link{distrib_penalty}} instantiates: a penalty built by
#' applying a univariate \pkg{distributions7} log-density coordinatewise to
#' \eqn{D\beta}.
#'
#' @inheritParams penalty
#' @param parent The \pkg{distributions7} object.
#' @param kinks The declared non-differentiable points of the parent's
#'   log-density in its argument.
#'
#' @return An object of class \code{DistribPenalty}.
#'
#' @seealso \code{\link{distrib_penalty}}
#' @examples
#' S7::S7_inherits(ridge_penalty(n_coef = 2), DistribPenalty)
#' @keywords internal
#' @export
DistribPenalty <- S7::new_class(
  name = "DistribPenalty",
  parent = penalty,
  properties = list(
    parent = S7::class_any,
    kinks = S7::class_numeric
  )
)

#' @title Construct a Separable Penalty From a Distribution
#'
#' @description
#' \eqn{\rho(D\beta;\theta) = -\sum_j \log f((D\beta)_j;\theta)} for a
#' univariate \pkg{distributions7} density \eqn{f}. The hyperparameters ARE
#' the distribution's free parameters; bounds and links are read off the
#' distribution object rather than restated, and every derivative is the
#' distribution's, reassembled: the gradient is \eqn{-D'\ell^{(y)}}, the
#' Hessian \eqn{-D'\mathrm{diag}(\ell^{(yy)})D}, the theta blocks the summed
#' score and Hessian, and the mixed block one
#' \eqn{-D'\ell^{(y\theta_k)}} per hyperparameter -- which is what
#' \code{\link[distributions7]{distrib_cross_y}} exists for.
#'
#' @details
#' Ridge, lasso and the heavy-tailed prior are this construction at a
#' \code{\link[distributions7]{fixed}} Gaussian, Laplace and Student t;
#' see \code{\link{ridge_penalty}}. The normalizing constant comes with the
#' density and is kept, so the value is exactly the negative log prior
#' density and a free scale or a free \eqn{\nu} is estimable.
#'
#' @param d A univariate continuous \pkg{distributions7} object; typically a
#'   \code{fixed()} wrapper holding the location at zero.
#' @param map The matrix \eqn{D}, or \code{NULL} (default) for the identity.
#' @param n_coef The number of coefficients; required when \code{map} is
#'   \code{NULL}, ignored otherwise.
#' @param kinks The points where the parent's log-density is not
#'   differentiable in its argument. \code{NULL}, the default, derives them
#'   from the parent with \code{\link{distrib_kinks}}; pass a numeric vector to
#'   say so directly, or \code{numeric(0)} to declare there are none.
#'
#' @return An object of class \code{DistribPenalty}.
#'
#' @examples
#' d <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
#' pen <- distrib_penalty(d, n_coef = 3)
#' penalty_value(pen, c(1, 0, -1), list(sigma = 2))
#'
#' @seealso \code{\link{quadratic_penalty}}, \code{\link{additive_penalty}}, \code{\link{structured_penalty}}
#' @export
distrib_penalty <- function(d, map = NULL, n_coef = NULL, kinks = NULL) {
  if (!S7::S7_inherits(d, distributions7::distrib)) {
    stop("'d' must be a distributions7 object.", call. = FALSE)
  }
  if (is.null(kinks)) kinks <- distrib_kinks(d)
  if (d@dimension != "univariate") {
    stop("'d' must be univariate.", call. = FALSE)
  }
  if (is.null(map)) {
    if (is.null(n_coef)) {
      stop("'n_coef' is required when 'map' is NULL.", call. = FALSE)
    }
    q <- as.integer(n_coef)
  } else {
  # A map that is already a Matrix is KEPT as it is: `as.matrix()` here would
  # densify a diagonal or sparse map, which is the whole cost the map exists
  # to avoid -- a diagonal one is a per-coordinate rescaling and costs q
  # numbers, its dense form q^2.
    map <- as_map(map)
    q <- ncol(map)
  }
  DistribPenalty(
    penalty_name = sprintf("separable [%s]", d@distrib_name),
    map = map,
    n_coef = q,
    params = d@params,
    params_bounds = d@params_bounds,
    link_params = d@link_params,
    params_smooth = stats::setNames(rep(TRUE, length(d@params)), d@params),
    parent = d,
    kinks = kinks
  )
}


#' Where the Parent's Log-Density Has a Kink
#'
#' @description
#' The points a separable penalty is not differentiable at, derived from the
#' distribution it is built on.
#'
#' @details
#' A distribution records which of its parameters the log-likelihood is
#' differentiable in, through \code{params_smooth}, and for a location family a
#' location that is not smooth is a kink in the argument at that location. A
#' penalty is the negative log-density read in the coefficient, so a
#' \code{\link[distributions7]{fixed}} wrapper holding such a parameter at a
#' value puts the kink there: \code{fixed(laplace_distrib(), mu = 0)} is the
#' lasso and has a kink at zero.
#'
#' Each candidate is the value its parameter is held at, and whether it is a
#' kink is then measured rather than inferred, by comparing the one-sided
#' derivatives of the log-density across it. Inferring alone would put a kink
#' on any family whose non-smooth parameter is not a location; measuring alone
#' would need somewhere to look. A candidate whose derivative does not jump is
#' dropped.
#'
#' A parent that declares every parameter smooth, or that fixes none of the
#' ones it declares non-smooth, has no candidate and gets \code{numeric(0)}.
#' Nothing is taken from a parameter that is free, its value being whatever the
#' hyperparameters say at the time.
#'
#' @param d A univariate \pkg{distributions7} object.
#'
#' @return A numeric vector, possibly empty.
#'
#' @seealso \code{\link{distrib_penalty}}, \code{\link{penalty_kinks}}
#'
#' @examples
#' distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
#'                                     mu = 0))
#' distrib_kinks(distributions7::fixed(distributions7::gaussian1_distrib(),
#'                                     mu = 0))
#'
#' @export
distrib_kinks <- function(d) {
  # the wrapper is recognized by what it carries rather than by its class:
  # the two fixed classes are not exported, and a family written outside
  # distributions7 may hold its parameters the same way
  if (!all(c("parent_distrib", "fixed_params") %in% S7::prop_names(d))) {
    return(numeric(0))
  }
  smooth <- d@parent_distrib@params_smooth
  fixed_at <- d@fixed_params
  rough <- names(smooth)[!smooth]
  cand <- unlist(fixed_at[intersect(rough, names(fixed_at))],
                 use.names = FALSE)
  cand <- unique(cand[is.finite(cand)])
  if (!length(cand)) return(numeric(0))
  th <- probe_theta(d)
  as.numeric(cand[vapply(cand, function(k) has_jump(d, th, k), logical(1))])
}


#' A Point to Read the Parent At
#'
#' @description
#' The midpoint of each free parameter's bounds, the probe rule the toolkit
#' already uses where a value is needed and none is in force.
#'
#' @param d A \pkg{distributions7} object.
#'
#' @return A named list.
#'
#' @keywords internal
probe_theta <- function(d) {
  stats::setNames(lapply(d@params, function(p) {
    b <- d@params_bounds[[p]]
    if (is.finite(b[1L]) && is.finite(b[2L])) (b[1L] + b[2L]) / 2
    else if (is.finite(b[1L])) b[1L] + 1
    else if (is.finite(b[2L])) b[2L] - 1
    else 1
  }), d@params)
}


#' Does the Log-Density's Slope Jump Across a Point?
#'
#' @description
#' Compares the derivative in the argument just either side of it.
#'
#' @param d A \pkg{distributions7} object.
#' @param theta Where to read it.
#' @param at The candidate point.
#' @param eps How far either side.
#'
#' @return A single logical.
#'
#' @keywords internal
has_jump <- function(d, theta, at, eps = 1e-5) {
  g <- tryCatch({
    up <- distributions7::distrib_grad_y(d, at + eps, theta)
    dn <- distributions7::distrib_grad_y(d, at - eps, theta)
    abs(as.numeric(up) - as.numeric(dn))
  }, error = function(e) NA_real_)
  isTRUE(is.finite(g) && g > 1e-6 * max(1, abs(at)))
}


#' @title Named Separable Penalties
#'
#' @description
#' The canonical instances of \code{\link{distrib_penalty}}, shipped as
#' constructors so the model layer can name what it means. Ridge is the
#' Gaussian at zero with the scale free; the lasso is the Laplace in location
#' and rate (`laplace2`) at zero, so the free hyperparameter is the rate
#' \eqn{\lambda} and the value is \eqn{\lambda\lVert D\beta\rVert_1} up to
#' its constant, with the kink declared; the heavy-tailed prior is the
#' Student t at zero,
#' whose \eqn{\nu} is estimable exactly because the normalizing constant is
#' kept. The elastic net is the product of the Laplace and the Gaussian at
#' zero, normalized (\code{\link[distributions7]{enet_distrib}}), so its
#' hyperparameters are the overall rate \eqn{\lambda} and the mixing weight
#' \eqn{\alpha} and its value is
#' \eqn{\lambda\{\alpha\lVert D\beta\rVert_1 +
#' (1-\alpha)\lVert D\beta\rVert_2^2/2\}} up to a constant. That constant
#' depends on both hyperparameters, which is what makes them estimable by a
#' marginal criterion and what a penalty written as a formula would not
#' have.
#'
#' @param map The matrix \eqn{D}, or \code{NULL} (default) for the identity.
#' @param n_coef The number of coefficients when \code{map} is \code{NULL}.
#'
#' @return An object of class \code{DistribPenalty}.
#'
#' @examples
#' pen <- ridge_penalty(n_coef = 2)
#' penalty_gradient(pen, c(1, -1), list(sigma = 1))
#'
#' @references
#' Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased
#' estimation for nonorthogonal problems. \emph{Technometrics} 12, 55-67.
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' \emph{Journal of the Royal Statistical Society, Series B} 58, 267-288.
#'
#' Zou, H. and Hastie, T. (2005). Regularization and variable selection via
#' the elastic net. \emph{Journal of the Royal Statistical Society, Series B}
#' 67, 301-320.
#'
#' @seealso \code{\link{distrib_penalty}}, \code{\link{scad_penalty}}, \code{\link{quadratic_penalty}}
#' @export
ridge_penalty <- function(map = NULL, n_coef = 1L) {
  distrib_penalty(
    distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
    map = map, n_coef = n_coef
  )
}

#' @rdname ridge_penalty
#' @export
lasso_penalty <- function(map = NULL, n_coef = 1L) {
  distrib_penalty(
    distributions7::fixed(distributions7::laplace2_distrib(), mu = 0),
    map = map, n_coef = n_coef, kinks = 0
  )
}

#' @rdname ridge_penalty
#' @export
elasticnet_penalty <- function(map = NULL, n_coef = 1L) {
  distrib_penalty(
    distributions7::fixed(distributions7::enet_distrib(), mu = 0),
    map = map, n_coef = n_coef, kinks = 0
  )
}

#' @rdname ridge_penalty
#' @export
heavy_penalty <- function(map = NULL, n_coef = 1L) {
  distrib_penalty(
    distributions7::fixed(distributions7::student_t1_distrib(), mu = 0),
    map = map, n_coef = n_coef
  )
}

#' @title Separable Penalty Methods
#' @name penalty_value.DistribPenalty
#' @description
#' The parent distribution's quantities, reassembled through the map; see
#' \code{\link{distrib_penalty}}.
#' @param pen A \code{DistribPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of the parent's free parameters.
#' @param scale Handled by the generic.
#' @param ... Unused.
#' @return See the generic pages.
#' @keywords internal
S7::method(penalty_value, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  -sum(distributions7::distrib_pdf(pen@parent, t, theta, log = TRUE))
}

#' @rdname penalty_value.DistribPenalty
#' @name penalty_gradient.DistribPenalty
#' @keywords internal
S7::method(penalty_gradient, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  -map_back(pen, distributions7::distrib_grad_y(pen@parent, t, theta) + 0 * t)
}

#' @rdname penalty_value.DistribPenalty
#' @name penalty_hessian.DistribPenalty
#' @keywords internal
S7::method(penalty_hessian, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  -map_quad(pen, distributions7::distrib_hess_y(pen@parent, t, theta) + 0 * t)
}

#' @rdname penalty_value.DistribPenalty
#' @name penalty_grad_theta.DistribPenalty
#' @keywords internal
S7::method(penalty_grad_theta, DistribPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    g <- distributions7::distrib_gradient(pen@parent, t, theta)
    lapply(g, function(v) -sum(v))
  }

#' @rdname penalty_value.DistribPenalty
#' @name penalty_hess_theta.DistribPenalty
#' @keywords internal
S7::method(penalty_hess_theta, DistribPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    H <- distributions7::distrib_hessian(pen@parent, t, theta)
    # reorder onto the diagonals-first convention of ptheta_pairs
    prs <- ptheta_pairs(pen@params)
    hn <- names(H)
    stats::setNames(lapply(names(prs), function(nm) {
      pr <- prs[[nm]]
      alt <- paste0(pr[2], "_", pr[1])
      v <- if (nm %in% hn) H[[nm]] else H[[alt]]
      -sum(v)
    }), names(prs))
  }

#' @rdname penalty_value.DistribPenalty
#' @name penalty_cross.DistribPenalty
#' @keywords internal
S7::method(penalty_cross, DistribPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    cy <- distributions7::distrib_cross_y(pen@parent, t, theta)
    lapply(cy, function(v) -map_back(pen, v + 0 * t))
  }

#' @rdname penalty_value.DistribPenalty
#' @name penalty_kinks.DistribPenalty
#' @keywords internal
S7::method(penalty_kinks, DistribPenalty) <- function(pen, theta, ...) {
  pen@kinks
}

#' @rdname penalty_value.DistribPenalty
#' @name is_proper.DistribPenalty
#' @keywords internal
S7::method(is_proper, DistribPenalty) <- function(pen, ...) {
  is.null(pen@map) || nrow(pen@map) >= ncol(pen@map)
}
