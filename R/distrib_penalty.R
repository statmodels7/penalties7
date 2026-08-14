#' @include generics.R
NULL

#' @title S7 Class for the Separable Penalty
#'
#' @description
#' The class \code{\link{distrib_penalty}} instantiates: a penalty built by
#' applying a \pkg{distributions7} log-density to the successive blocks of
#' \eqn{D\beta}, the block being one coordinate for a univariate parent and
#' \eqn{p} of them for a \eqn{p}-variate one.
#'
#' @inheritParams penalty
#' @param parent The \pkg{distributions7} object.
#' @param kinks The declared non-differentiable points of the parent's
#'   log-density in its argument.
#' @param block The block width, \code{1} for a univariate parent and the
#'   parent's dimension otherwise.
#'
#' @return An object of class \code{DistribPenalty}.
#'
#' @seealso \code{\link{distrib_penalty}}
#' @examples
#' S7::S7_inherits(lasso_penalty(n_coef = 2), DistribPenalty)
#' @keywords internal
#' @export
DistribPenalty <- S7::new_class(
  name = "DistribPenalty",
  parent = penalty,
  properties = list(
    parent = S7::class_any,
    kinks = S7::class_numeric,
    block = S7::class_integer
  )
)

#' @title Construct a Separable Penalty From a Distribution
#'
#' @description
#' \eqn{\rho(D\beta;\theta) = -\sum_i \log f(b_i;\theta)} for a
#' \pkg{distributions7} density \eqn{f} read at the successive blocks
#' \eqn{b_i} of \eqn{D\beta}. A univariate parent gives blocks of one
#' coordinate, which is the separable penalty; a \eqn{p}-variate parent gives
#' blocks of \eqn{p}, which is a prior that lets the coordinates of one block
#' depend on each other while the blocks stay independent.
#'
#' The hyperparameters ARE
#' the distribution's free parameters; bounds and links are read off the
#' distribution object rather than restated, and every derivative is the
#' distribution's, reassembled: the gradient is \eqn{-D'\ell^{(y)}}, the
#' Hessian \eqn{-D'\mathrm{diag}(\ell^{(yy)})D} (block diagonal rather than
#' diagonal when the parent is multivariate), the theta blocks the summed
#' score and Hessian, and the mixed block one
#' \eqn{-D'\ell^{(y\theta_k)}} per hyperparameter -- which is what
#' \code{\link[distributions7]{distrib_cross_y}} exists for.
#'
#' @details
#' Lasso, the elastic net and the heavy-tailed prior are this construction at
#' a \code{\link[distributions7]{fixed}} Laplace, elastic net and Student t;
#' see \code{\link{ridge_penalty}}. The normalizing constant comes with the
#' density and is kept, so the value is exactly the negative log prior
#' density and a free scale or a free \eqn{\nu} is estimable.
#'
#' A multivariate parent is centered by the caller, typically through
#' \code{\link[distributions7]{fixed}} at a zero mean, and its matrix
#' parameter carries the dependence within a block. It has no proximal
#' operator: that operator acts one coordinate at a time and the coordinates
#' of a block do not separate.
#'
#' @param d A continuous \pkg{distributions7} object; typically a
#'   \code{fixed()} wrapper holding the location at zero. A multivariate
#'   parent of dimension \eqn{p} is read blockwise, and the number of
#'   coefficients must then be a multiple of \eqn{p}.
#' @param map The matrix \eqn{D}, or \code{NULL} (default) for the identity.
#' @param n_coef The number of coefficients; required when \code{map} is
#'   \code{NULL}, ignored otherwise.
#' @param kinks The points where the parent's log-density is not
#'   differentiable in its argument. \code{NULL}, the default, derives them
#'   from the parent with \code{\link{distrib_kinks}}; pass a numeric vector to
#'   say so directly, or \code{numeric(0)} to declare there are none. A
#'   multivariate parent has none: a kink is a point of a scalar argument.
#'
#' @return An object of class \code{DistribPenalty}.
#'
#' @examples
#' d <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
#' pen <- distrib_penalty(d, n_coef = 3)
#' penalty_value(pen, c(1, 0, -1), list(sigma = 2))
#'
#' # a correlated prior over two coefficients per block, three blocks
#' mv <- distributions7::fixed(distributions7::mvgaussian_distrib(2),
#'                             mu1 = 0, mu2 = 0)
#' pen2 <- distrib_penalty(mv, n_coef = 6)
#' penalty_value(pen2, c(1, 0, -1, 0.5, 0.2, -0.3),
#'               list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.4))
#'
#' @seealso \code{\link{quadratic_penalty}}, \code{\link{additive_penalty}}, \code{\link{structured_penalty}}
#' @export
distrib_penalty <- function(d, map = NULL, n_coef = NULL, kinks = NULL) {
  if (!S7::S7_inherits(d, distributions7::distrib)) {
    stop("'d' must be a distributions7 object.", call. = FALSE)
  }
  mv <- identical(d@dimension, "multivariate")
  if (!mv && !identical(d@dimension, "univariate")) {
    stop(sprintf("'d' is %s, and a penalty is read one block at a time.",
                 d@dimension), call. = FALSE)
  }
  # A kink is a point of a scalar argument, so a multivariate parent has no
  # candidate: this is said here rather than inherited from the univariate
  # route by accident.
  if (is.null(kinks)) kinks <- if (mv) numeric(0) else distrib_kinks(d)
  blk <- if (mv) as.integer(d@n_dim) else 1L
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
  rows <- if (is.null(map)) q else nrow(map)
  if (rows %% blk != 0L) {
    stop(sprintf(paste0(
      "'%s' is %d-variate, so the %d values it is read at must divide into\n",
      "  whole blocks; %d does not."),
      d@distrib_name, blk, rows, rows), call. = FALSE)
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
    kinks = kinks,
    block = blk
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
#' constructors so the model layer can name what it means. Each is written
#' on the chart whose hyperparameter MEASURES THE SHRINKAGE, so that a larger
#' value shrinks harder in all of them. Ridge is the exception to the branch
#' rather than to the rule: it is the Gaussian prior at zero, and that prior
#' written by its PRECISION is exactly the quadratic penalty at the identity,
#' the same value to the last bit, so it is built there and its
#' hyperparameter is the lambda that branch already carries -- one name for
#' one number. The lasso is the Laplace in location
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
#' penalty_gradient(pen, c(1, -1), list(lambda = 1))
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
  k <- if (is.null(map)) as.integer(n_coef) else nrow(as.matrix(map))
  quadratic_penalty(diag(1, k), map = map)
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
NULL

#' The Parent's Argument, Shaped and Unshaped
#'
#' @description
#' \code{dp_arg} reshapes \eqn{D\beta} into the argument the parent reads --
#' the vector itself for a univariate parent, one row per block for a
#' multivariate one -- and \code{dp_flat} undoes it.
#'
#' @details
#' The blocks are the SUCCESSIVE stretches of \eqn{D\beta}, so the reshaping
#' fills by row: block \eqn{i} occupies positions \eqn{(i-1)p+1, \dots, ip}.
#' That is the order a grouped design assembles its coefficients in, and the
#' order \eqn{I_m \otimes \Sigma} assumes.
#'
#' @param pen A \code{DistribPenalty} object.
#' @param t The mapped coefficient vector.
#' @param g The parent's answer, a vector or a matrix of one row per block.
#'
#' @return A vector or a matrix.
#'
#' @keywords internal
dp_arg <- function(pen, t) {
  if (pen@block == 1L) t else matrix(t, ncol = pen@block, byrow = TRUE)
}

#' @rdname dp_arg
#' @keywords internal
dp_flat <- function(pen, g) {
  if (pen@block == 1L) as.numeric(g) else as.numeric(t(g))
}

#' The Block-Diagonal Middle Matrix of a Blockwise Parent
#'
#' @description
#' Assembles \eqn{\partial^2\ell/\partial b\partial b'} over the blocks. The
#' parent returns one \eqn{p \times p} matrix when its response Hessian does
#' not depend on the observation, as the gaussian's does not, and a
#' \eqn{p \times p \times n} array when it does.
#'
#' @param pen A \code{DistribPenalty} object.
#' @param h The parent's \code{distrib_hess_y}.
#' @param nblk The number of blocks.
#'
#' @return A symmetric matrix of side \code{nblk * pen@block}.
#'
#' @keywords internal
dp_blockdiag <- function(pen, h, nblk) {
  p <- pen@block
  out <- matrix(0, nblk * p, nblk * p)
  const <- is.matrix(h)
  for (i in seq_len(nblk)) {
    ix <- (i - 1L) * p + seq_len(p)
    out[ix, ix] <- if (const) h else h[, , i]
  }
  out
}

S7::method(penalty_value, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  -sum(distributions7::distrib_pdf(pen@parent, dp_arg(pen, t), theta,
                                   log = TRUE))
}

#' @rdname penalty_value.DistribPenalty
#' @name penalty_gradient.DistribPenalty
#' @keywords internal
S7::method(penalty_gradient, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  a <- dp_arg(pen, t)
  g <- distributions7::distrib_grad_y(pen@parent, a, theta) + 0 * a
  -map_back(pen, dp_flat(pen, g))
}

#' @rdname penalty_value.DistribPenalty
#' @name penalty_hessian.DistribPenalty
#' @keywords internal
S7::method(penalty_hessian, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  if (pen@block == 1L) {
    return(-map_quad(pen,
      distributions7::distrib_hess_y(pen@parent, t, theta) + 0 * t))
  }
  a <- dp_arg(pen, t)
  h <- distributions7::distrib_hess_y(pen@parent, a, theta)
  -map_quad_full(pen, dp_blockdiag(pen, h, nrow(a)))
}

#' @rdname penalty_value.DistribPenalty
#' @name penalty_grad_theta.DistribPenalty
#' @keywords internal
S7::method(penalty_grad_theta, DistribPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    g <- distributions7::distrib_gradient(pen@parent, dp_arg(pen, t), theta)
    lapply(g, function(v) -sum(v))
  }

#' @rdname penalty_value.DistribPenalty
#' @name penalty_hess_theta.DistribPenalty
#' @keywords internal
S7::method(penalty_hess_theta, DistribPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    H <- distributions7::distrib_hessian(pen@parent, dp_arg(pen, t), theta)
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
    a <- dp_arg(pen, t)
    cy <- distributions7::distrib_cross_y(pen@parent, a, theta)
    lapply(cy, function(v) -map_back(pen, dp_flat(pen, v + 0 * a)))
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

#' @rdname penalty_value.DistribPenalty
#' @name penalty_readable.DistribPenalty
#' @keywords internal
S7::method(penalty_readable, DistribPenalty) <- function(pen, theta, ...) {
  # only a blockwise parent has coordinates that are not the quantities: a
  # univariate one's hyperparameters are a scale, a rate, a shape, each read
  # on its own scale already
  if (pen@block == 1L) return(NULL)
  distributions7::mv_derived(pen@parent, theta)
}
