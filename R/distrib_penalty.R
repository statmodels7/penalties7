#' @include generics.R
NULL

#' @title S7 Class for the Separable Penalty
#'
#' @description
#' The class [distrib_penalty()] builds: a penalty formed by applying a
#' \pkg{distributions7} log-density to the successive blocks of \eqn{D\beta}.
#' Beyond the seven properties every penalty carries it stores the parent
#' distribution, the points where the parent's log-density has a kink, and the
#' block width, which is one for a univariate parent and the dimension for a
#' multivariate one.
#'
#' @details
#' The hyperparameters are the parent's own free parameters, and `params`,
#' `params_bounds` and `link_params` are read straight off the distribution
#' object. A caller who knows the parent therefore knows the penalty's
#' hyperparameters, their bounds and their links without being told again.
#'
#' @inheritParams penalty
#' @param parent The \pkg{distributions7} object the log-density comes from,
#'   usually a [distributions7::fixed()] wrapper holding the location at zero.
#' @param kinks A numeric vector of the points where the parent's log-density
#'   is not differentiable in its argument, possibly empty. Derived by
#'   [distrib_kinks()] unless the constructor was told otherwise.
#' @param block The block width: `1L` for a univariate parent and the parent's
#'   `n_dim` for a multivariate one.
#'
#' @return An S7 object of class `DistribPenalty`, inheriting from [penalty()],
#'   with the seven inherited properties and the three above.
#'
#' @seealso [distrib_penalty()] for the constructor to use,
#'   [ridge_penalty()] and its siblings for the named instances,
#'   [penalty_value.DistribPenalty()] for what the branch computes.
#'
#' @examples
#' pen <- lasso_penalty(n_coef = 2)
#' S7::S7_inherits(pen, DistribPenalty)
#'
#' # The parent is a fixed Laplace, the kink is at zero, blocks are single
#' # coordinates, and the hyperparameter is the parent's own rate.
#' pen@parent
#' pen@kinks
#' pen@block
#' pen@params
#'
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
#' Builds the penalty \eqn{\rho(D\beta;\theta) = -\sum_i \log f(b_i;\theta)},
#' where \eqn{f} is a \pkg{distributions7} density and the \eqn{b_i} are the
#' successive blocks of \eqn{D\beta}. A univariate parent gives blocks of one
#' coordinate, which is the separable penalty; a \eqn{p}-variate parent gives
#' blocks of \eqn{p}, a prior under which the coordinates of one block depend
#' on each other while the blocks stay independent.
#'
#' @details
#' # The hyperparameters are the distribution's
#'
#' `params`, `params_bounds` and `link_params` are read off the distribution
#' rather than restated, and every derivative is the distribution's own,
#' reassembled through the map:
#'
#' \deqn{\frac{\partial\rho}{\partial\beta} = -D'\ell^{(y)}, \qquad
#'   \frac{\partial^2\rho}{\partial\beta^2}
#'     = -D'\,\mathrm{diag}(\ell^{(yy)})\,D, \qquad
#'   \frac{\partial^2\rho}{\partial\beta\,\partial\theta_k}
#'     = -D'\ell^{(y\theta_k)},}
#'
#' with the middle matrix block diagonal instead of diagonal when the parent is
#' multivariate. The hyperparameter blocks are the parent's score and Hessian
#' summed over the blocks, and the mixed block is
#' [distributions7::distrib_cross_y()], which exists for this.
#'
#' # The constant is kept
#'
#' The normalizing constant comes with the density and is not dropped, so the
#' value is exactly the negative log prior density. A free scale or a free
#' \eqn{\nu} is then estimable: with the constant dropped, a prior scale could
#' be sent to infinity for nothing.
#'
#' # A multivariate parent
#'
#' Centering is the caller's, typically through [distributions7::fixed()] at a
#' zero mean, and the matrix parameter carries the dependence within a block.
#' The number of values the penalty is read at must divide into whole blocks,
#' and the constructor rejects a count that does not, naming both numbers.
#'
#' A blockwise penalty has **no proximal operator**: that operator acts one
#' coordinate at a time and the coordinates of a block do not separate.
#' [has_prox()] answers `FALSE` and [penalty_prox()] rejects.
#'
#' @param d A continuous \pkg{distributions7} object, univariate or
#'   multivariate; typically a [distributions7::fixed()] wrapper holding the
#'   location at zero. A discrete or otherwise non-continuous object is
#'   rejected, naming its `dimension`.
#' @param map The matrix \eqn{D}, with one column per coefficient, or `NULL`
#'   (the default) for the identity. A \pkg{Matrix} object keeps its own
#'   storage, a diagonal map being what standardization is.
#' @param n_coef The number of coefficients. Required when `map` is `NULL` and
#'   ignored otherwise, `ncol(map)` being the count then.
#' @param kinks The points where the parent's log-density is not
#'   differentiable in its argument. `NULL`, the default, derives them with
#'   [distrib_kinks()]; pass a numeric vector to say so directly, or
#'   `numeric(0)` to declare there are none. A multivariate parent is given
#'   `numeric(0)` without asking, a kink being a point of a scalar argument.
#'
#' @return A [DistribPenalty()] object whose hyperparameters, bounds and links
#'   are the parent distribution's.
#'
#' @examples
#' # A Gaussian prior at zero: the separable twin of the ridge, written by
#' # its standard deviation rather than by a precision.
#' d <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
#' pen <- distrib_penalty(d, n_coef = 3)
#' pen@params
#' penalty_value(pen, c(1, 0, -1), list(sigma = 2))
#'
#' # Which is exactly the negative log-density of that prior.
#' -sum(stats::dnorm(c(1, 0, -1), sd = 2, log = TRUE))
#'
#' # A correlated prior over two coefficients per block, three blocks. The
#' # hyperparameters are the matrix parameter's free values.
#' mv <- distributions7::fixed(distributions7::mvgaussian1_distrib(2),
#'                             mu1 = 0, mu2 = 0)
#' pen2 <- distrib_penalty(mv, n_coef = 6)
#' pen2@block
#' pen2@params
#' penalty_value(pen2, c(1, 0, -1, 0.5, 0.2, -0.3),
#'               list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.4))
#'
#' # And it has no proximal operator, the coordinates of a block being tied.
#' has_prox(pen2)
#'
#' # A count that does not divide into blocks is rejected.
#' try(distrib_penalty(mv, n_coef = 5))
#'
#' @seealso [ridge_penalty()], [lasso_penalty()], [elasticnet_penalty()] and
#'   [heavy_penalty()] for the named instances,
#'   [quadratic_penalty()] and [structured_penalty()] for the quadratic
#'   branches, [distrib_kinks()] for how the kinks are found,
#'   [distributions7::fixed()] for the wrapper that centers a parent.
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
#' Returns the points at which a separable penalty built on this distribution
#' would not be differentiable, derived from the distribution itself. Returns
#' `numeric(0)` where there are none, which is the answer for every family the
#' toolkit ships except the Laplace and the elastic net.
#'
#' @details
#' # How a candidate is found
#'
#' A distribution records which of its parameters the log-likelihood is
#' differentiable in, through `params_smooth`, and for a location family a
#' location that is not smooth is a kink in the argument at that location. A
#' penalty is the negative log-density read in the coefficient, so a
#' [distributions7::fixed()] wrapper holding such a parameter at a value puts
#' the kink there: `fixed(laplace_distrib(), mu = 0)` is the lasso and has a
#' kink at zero, while `fixed(laplace_distrib(), mu = 2)` has one at two.
#'
#' Nothing is taken from a parameter that is **free**, its value being whatever
#' the hyperparameters say at the time. A Laplace with an unfixed location
#' therefore reports no kink at all.
#'
#' # And how it is confirmed
#'
#' Each candidate is the value its parameter is held at, and whether it is a
#' kink is then measured rather than inferred, by comparing the one-sided
#' derivatives of the log-density across it at a probe point
#' ([probe_theta()], [has_jump()]). Inferring alone would put a kink on any
#' family whose non-smooth parameter is not a location; measuring alone would
#' need somewhere to look. A candidate whose derivative does not jump is
#' dropped.
#'
#' The wrapper is recognized by the properties it carries, `parent_distrib`
#' and `fixed_params`, and not by its class: the two `fixed` classes are not
#' exported, and a family written outside \pkg{distributions7} may hold its
#' parameters the same way.
#'
#' @param d A univariate \pkg{distributions7} object. An object carrying
#'   neither `parent_distrib` nor `fixed_params` gets `numeric(0)` without
#'   further inspection.
#'
#' @return A numeric vector of the confirmed kink points, possibly empty, with
#'   duplicates and non-finite candidates removed.
#'
#' @seealso [distrib_penalty()], which calls this unless told the kinks
#'   directly, [penalty_kinks()] for what a built penalty reports,
#'   [has_jump()] for the measurement.
#'
#' @examples
#' # The lasso's parent: a kink at the location it is held at.
#' distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
#'                                     mu = 0))
#' distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
#'                                     mu = 2))
#'
#' # A Gaussian is smooth everywhere, so there is no candidate to confirm.
#' distrib_kinks(distributions7::fixed(distributions7::gaussian1_distrib(),
#'                                     mu = 0))
#'
#' # And a free location gives nothing, its value not being settled.
#' distrib_kinks(distributions7::laplace_distrib())
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
#' Returns one value per free parameter, at which the parent can be evaluated
#' when no hyperparameters are in force. Used by [distrib_kinks()], which has
#' to differentiate the log-density before any penalty exists.
#'
#' @details
#' The rule is the toolkit's probe rule: the midpoint of a parameter's bounds
#' where both are finite, one above the lower bound where only that is finite,
#' one below the upper where only that is, and `1` where neither is. So a scale
#' on \eqn{(0, \infty)} is probed at 1 and a mixing weight on \eqn{(0, 1)} at
#' 0.5.
#'
#' Whether a point is a kink does not depend on where the other parameters sit,
#' so any admissible point serves and the rule only has to give one.
#'
#' @param d A \pkg{distributions7} object.
#'
#' @return A named list of one number per element of `d@params`, in that order.
#'   `list()` for a distribution with no free parameters.
#'
#' @seealso [distrib_kinks()], [has_jump()]
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
#' Compares the derivative of the log-density in its argument just either side
#' of a candidate point, and answers `TRUE` when the two differ by more than a
#' relative threshold. This is the measurement that turns a candidate kink,
#' inferred from `params_smooth`, into a confirmed one.
#'
#' @details
#' The comparison is `abs(grad_y(at + eps) - grad_y(at - eps))` against
#' `1e-6 * max(1, abs(at))`, so the threshold is absolute near the origin and
#' relative away from it. At a genuine kink the difference is of the order of
#' the jump itself and does not shrink with `eps`: for a Laplace at zero it is
#' twice the rate. At a smooth point it is \eqn{O(\varepsilon)} and falls well
#' below the threshold.
#'
#' A parent that cannot be differentiated at either point returns `NA` from
#' [distributions7::distrib_grad_y()] or raises, and both are caught and read
#' as `FALSE`: an unconfirmed candidate is dropped rather than assumed.
#'
#' @param d A \pkg{distributions7} object.
#' @param theta A named list of parameter values to read it at, as
#'   [probe_theta()] returns.
#' @param at The candidate point, a single number.
#' @param eps How far either side to look, `1e-5` by default. Large enough that
#'   the one-sided derivatives are not dominated by rounding, small enough that
#'   a second feature does not fall inside the window.
#'
#' @return A single logical, `TRUE` when the slope jumps.
#'
#' @seealso [distrib_kinks()], [probe_theta()]
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


#' @title The Named Penalties
#'
#' @description
#' The four penalties a modeling layer names directly: `ridge_penalty()`,
#' `lasso_penalty()`, `elasticnet_penalty()` and `heavy_penalty()`. Each is a
#' particular prior centered at zero, and each is written on the chart whose
#' hyperparameter **measures the shrinkage**, so that a larger value shrinks
#' harder in all four.
#'
#' @details
#' # Ridge, which is on the other branch
#'
#' `ridge_penalty()` is the Gaussian prior at zero, and that prior written by
#' its **precision** is exactly [quadratic_penalty()] at the identity matrix,
#' the same value to the last bit. It is therefore built there, and its
#' hyperparameter is the `lambda` that branch already carries: one name for one
#' number. It is the only one of the four that is not a `DistribPenalty`, and
#' [is_quadratic()] answers `TRUE` for it alone.
#'
#' The separable twin still exists and is one line,
#' `distrib_penalty(fixed(gaussian1_distrib(), mu = 0), n_coef = q)`, whose
#' hyperparameter is the standard deviation \eqn{\sigma} instead. Reach for it
#' when a standard deviation is the quantity you want reported.
#'
#' # Lasso
#'
#' The Laplace in location and rate (`laplace2`) held at zero, so the free
#' hyperparameter is the rate \eqn{\lambda} and
#'
#' \deqn{\rho(\beta;\lambda) = \lambda\lVert D\beta\rVert_1
#'   - q\log(\lambda/2),}
#'
#' with \eqn{q} the number of values the penalty is read at. The kink at zero
#' is declared, so [penalty_kinks()] reports it and [penalty_prox()] gives the
#' soft threshold.
#'
#' # Elastic net
#'
#' The product of the Laplace and the Gaussian at zero, normalized
#' ([distributions7::enet_distrib()]), so the hyperparameters are the overall
#' rate \eqn{\lambda} and the mixing weight \eqn{\alpha} on \eqn{(0, 1)}:
#'
#' \deqn{\rho(\beta;\lambda,\alpha) = \lambda\left\{
#'   \alpha\lVert D\beta\rVert_1
#'   + (1-\alpha)\lVert D\beta\rVert_2^2/2\right\} + c(\lambda, \alpha).}
#'
#' The constant depends on **both** hyperparameters, which is why a marginal
#' criterion can estimate them: measured at
#' \eqn{(\lambda, \alpha)} of \eqn{(1, 0.5)}, \eqn{(2, 0.5)} and
#' \eqn{(1, 0.9)}, the per-coordinate constant is `0.7805`, `0.2711` and
#' `0.7001`. A penalty written as a formula, with no constant, would have
#' neither hyperparameter identified.
#'
#' # The heavy-tailed prior
#'
#' The Student t at zero, with a free scale \eqn{\sigma} and a free
#' \eqn{\nu}. Its \eqn{\nu} is estimable exactly because the normalizing
#' constant is kept, and estimating it is the point: a small \eqn{\nu} shrinks
#' small coefficients like a Gaussian prior and leaves large ones nearly alone,
#' which a Gaussian cannot do at any scale. Alone among the four it declares no
#' kink, so it produces no exact zeros.
#'
#' @param map The matrix \eqn{D}, with one column per coefficient, or `NULL`
#'   (the default) for the identity. Given a map, the number of coefficients is
#'   `ncol(map)` and `n_coef` is ignored.
#' @param n_coef The number of coefficients when `map` is `NULL`. A single
#'   whole number, `1L` by default.
#'
#' @return `ridge_penalty()` a [QuadraticPenalty()] with the hyperparameter
#'   `lambda`, a precision, on \eqn{(0, \infty)}.
#'   `lasso_penalty()` a [DistribPenalty()] with `lambda`, a rate, on
#'   \eqn{(0, \infty)} and a kink at zero.
#'   `elasticnet_penalty()` a [DistribPenalty()] with `lambda` on
#'   \eqn{(0, \infty)} and `alpha` on \eqn{(0, 1)}, and a kink at zero.
#'   `heavy_penalty()` a [DistribPenalty()] with `sigma` and `nu` both on
#'   \eqn{(0, \infty)}, and no kink.
#'
#' @references
#' Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased estimation
#' for nonorthogonal problems. *Technometrics* **12**, 55-67.
#'
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' *Journal of the Royal Statistical Society, Series B* **58**, 267-288.
#'
#' Zou, H. and Hastie, T. (2005). Regularization and variable selection via
#' the elastic net. *Journal of the Royal Statistical Society, Series B*
#' **67**, 301-320.
#'
#' @examples
#' b <- c(1, -0.5, 0.3)
#'
#' # Ridge is on the quadratic branch, so its lambda is a precision.
#' r <- ridge_penalty(n_coef = 3)
#' class(r)[1]
#' is_quadratic(r)
#' penalty_gradient(r, b, list(lambda = 2))
#'
#' # The lasso's value is lambda times the L1 norm, plus its constant.
#' l <- lasso_penalty(n_coef = 3)
#' penalty_value(l, b, list(lambda = 0.5))
#' 0.5 * sum(abs(b)) - 3 * log(0.5 / 2)
#'
#' # Only the lasso and the elastic net have a kink, so only they can set a
#' # coefficient exactly to zero.
#' penalty_kinks(l, list(lambda = 1))
#' penalty_kinks(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#'
#' # The heavy-tailed prior shrinks a large coefficient far less than a
#' # ridge at a comparable scale, and the small ones about as much.
#' penalty_prox(heavy_penalty(n_coef = 3), c(4, 0.3, -1), 1,
#'              list(sigma = 1, nu = 3))
#' penalty_prox(r, c(4, 0.3, -1), 1, list(lambda = 1))
#'
#' @seealso [distrib_penalty()] for the construction behind three of them,
#'   [quadratic_penalty()] for the one behind ridge,
#'   [scad_penalty()] and [mcp_penalty()] for the non-convex alternatives,
#'   [penalty_prox()] for the operators these have.
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

#' The Parent's Argument, Shaped and Unshaped
#'
#' @description
#' `dp_arg()` reshapes \eqn{D\beta} into the argument the parent reads: the
#' vector itself for a univariate parent, and a matrix of one row per block for
#' a multivariate one. `dp_flat()` undoes the reshaping, carrying the parent's
#' answer back to a vector in coefficient order.
#'
#' @details
#' The blocks are the **successive** stretches of \eqn{D\beta}, so the
#' reshaping fills by row: block \eqn{i} occupies positions
#' \eqn{(i-1)p+1, \dots, ip}. That is the order a grouped design assembles its
#' coefficients in, and the order \eqn{I_m \otimes \Sigma} assumes, so the two
#' agree without either being told about the other.
#'
#' Both are the identity when the block width is one, which is every univariate
#' parent, and the pair costs nothing there.
#'
#' @param pen A [DistribPenalty()] object.
#' @param t The mapped coefficient vector \eqn{D\beta}, of length a multiple of
#'   `pen@block`. `dp_arg()` only.
#' @param g The parent's answer: a vector for a univariate parent, a matrix of
#'   one row per block for a multivariate one. `dp_flat()` only.
#'
#' @return `dp_arg()` the vector unchanged when `pen@block` is `1L`, and
#'   otherwise a matrix with `pen@block` columns filled by row.
#'   `dp_flat()` a numeric vector, in coefficient order.
#'
#' @seealso [distrib_penalty()], [dp_blockdiag()] for the matrix-valued
#'   counterpart.
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
#' Assembles the parent's \eqn{\partial^2\ell/\partial b\partial b'} into the
#' block-diagonal matrix a multivariate separable penalty's Hessian is
#' sandwiched with: one \eqn{p \times p} block per block of coefficients, and
#' zeros between them.
#'
#' @details
#' The zeros between the blocks are the separability. Coordinates in different
#' blocks are independent under the prior, so the Hessian has no entry linking
#' them, however dependent the coordinates within a block are.
#'
#' Two input shapes are accepted because \pkg{distributions7} returns two. A
#' family whose response Hessian does not depend on the observation, as the
#' multivariate gaussian's does not, returns one \eqn{p \times p} matrix and it
#' is placed in every block; a family whose does returns a
#' \eqn{p \times p \times n} array and each slice goes to its own block. The
#' two are told apart with `is.matrix()`.
#'
#' @param pen A [DistribPenalty()] object with `pen@block` above one.
#' @param h The parent's [distributions7::distrib_hess_y()]: a
#'   \eqn{p \times p} matrix or a \eqn{p \times p \times n} array.
#' @param nblk The number of blocks, a single whole number.
#'
#' @return A symmetric base matrix of side `nblk * pen@block`, block diagonal
#'   with `nblk` blocks.
#'
#' @seealso [dp_arg()], and [map_quad_full()], which carries this through the
#'   map.
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

#' @title Value of a Separable Penalty
#' @name penalty_value.DistribPenalty
#'
#' @description
#' Returns \eqn{-\sum_i \log f(b_i;\theta)}, the negative log-density of the
#' parent summed over the blocks of \eqn{D\beta}. Because the parent supplies
#' its own normalizing constant, the value is exactly the negative log prior
#' density and needs nothing added.
#'
#' @details
#' The map is applied, the result reshaped into the parent's argument by
#' [dp_arg()], and one call to [distributions7::distrib_pdf()] with
#' `log = TRUE` gives every block's contribution at once. A univariate parent
#' sees the vector; a \eqn{p}-variate one sees a matrix of one row per block.
#'
#' @param pen A [DistribPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A named list of the parent's free parameters, already aligned
#'   and bound-checked by the generic against the parent's own bounds.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number.
#'
#' @examples
#' # A Gaussian prior at zero, which is the negative log-density of one.
#' pen <- distrib_penalty(
#'   distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
#'   n_coef = 3)
#' b <- c(1, 0, -1)
#' penalty_value(pen, b, list(sigma = 2))
#' -sum(stats::dnorm(b, sd = 2, log = TRUE))
#'
#' # The lasso's value is lambda times the L1 norm plus -q log(lambda/2).
#' penalty_value(lasso_penalty(n_coef = 3), b, list(lambda = 0.5))
#' 0.5 * sum(abs(b)) - 3 * log(0.5 / 2)
#'
#' @seealso [distrib_penalty()] for the construction,
#'   [penalty_gradient.DistribPenalty()] for the coefficient derivatives,
#'   [ridge_penalty()] and its siblings for the named instances.
#' @keywords internal
S7::method(penalty_value, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  -sum(distributions7::distrib_pdf(pen@parent, dp_arg(pen, t), theta,
                                   log = TRUE))
}

#' @title Coefficient Derivatives of a Separable Penalty
#' @name penalty_gradient.DistribPenalty
#'
#' @description
#' The parent's response derivatives, negated and carried back through the map.
#' `penalty_gradient()` returns \eqn{-D'\ell^{(y)}} and `penalty_hessian()`
#' returns \eqn{-D'\,\mathrm{diag}(\ell^{(yy)})\,D}, with the middle matrix
#' block diagonal instead of diagonal when the parent is multivariate.
#'
#' @details
#' Everything comes from [distributions7::distrib_grad_y()] and
#' [distributions7::distrib_hess_y()], which are closed form for every
#' continuous family, so no derivative is taken here.
#'
#' The `+ 0 * a` in both bodies is a recycling guard. A family whose derivative
#' does not depend on the argument, a Gaussian's second one for instance,
#' returns a single number where a vector is wanted, and adding zero times the
#' argument widens it without changing a value.
#'
#' At a kink the gradient is one element of the subdifferential, chosen by the
#' parent. For the lasso at \eqn{\beta_j = 0} that is `0`, where the true
#' subdifferential is \eqn{[-\lambda, \lambda]}; [penalty_prox()] is what
#' produces the exact zero.
#'
#' @param pen A [DistribPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of the parent's free parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` a numeric vector of length `pen@n_coef`;
#'   `penalty_hessian()` a symmetric base matrix of that side, diagonal under
#'   the identity map with a univariate parent and block diagonal with a
#'   multivariate one.
#'
#' @examples
#' b <- c(1, 0, -1)
#'
#' # A Gaussian prior at zero: the gradient is b/sigma^2 and the Hessian is
#' # constant.
#' pen <- distrib_penalty(
#'   distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
#'   n_coef = 3)
#' penalty_gradient(pen, b, list(sigma = 2))
#' b / 4
#' diag(penalty_hessian(pen, b, list(sigma = 2)))
#'
#' # The lasso's gradient is lambda times a sign, and 0 at the kink.
#' penalty_gradient(lasso_penalty(n_coef = 3), b, list(lambda = 1.5))
#'
#' # A multivariate parent gives a block-diagonal Hessian: dependence within
#' # a block, none between blocks.
#' mv <- distrib_penalty(
#'   distributions7::fixed(distributions7::mvgaussian1_distrib(2),
#'                         mu1 = 0, mu2 = 0), n_coef = 4)
#' round(penalty_hessian(mv, c(1, 0, -1, 0.5),
#'                       list(sigma_log_L1 = 0, sigma_log_L2 = 0,
#'                            sigma_L2.1 = 0.4)), 3)
#'
#' @seealso [penalty_value.DistribPenalty()] for the quantity differentiated,
#'   [penalty_grad_theta.DistribPenalty()] for the hyperparameter blocks,
#'   [penalty_prox()] for the step that produces exact zeros,
#'   [dp_blockdiag()] for the multivariate middle matrix.
#' @keywords internal
S7::method(penalty_gradient, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  a <- dp_arg(pen, t)
  g <- distributions7::distrib_grad_y(pen@parent, a, theta) + 0 * a
  -map_back(pen, dp_flat(pen, g))
}

#' @rdname penalty_gradient.DistribPenalty
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

#' @title Hyperparameter Derivatives of a Separable Penalty
#' @name penalty_grad_theta.DistribPenalty
#'
#' @description
#' The parent's own score and information in its parameters, negated and summed
#' over the blocks, plus the mixed block.
#' `penalty_grad_theta()` returns \eqn{-\sum_i \ell^{(\theta_k)}},
#' `penalty_hess_theta()` returns \eqn{-\sum_i \ell^{(\theta_k\theta_l)}}, and
#' `penalty_cross()` returns \eqn{-D'\ell^{(y\theta_k)}}, one coefficient
#' vector per hyperparameter.
#'
#' @details
#' The first two are [distributions7::distrib_gradient()] and
#' [distributions7::distrib_hessian()] read at the blocks and summed, the map
#' not entering: a hyperparameter of the prior does not act through \eqn{D}.
#' The third is [distributions7::distrib_cross_y()], which exists for this,
#' carried back through the map like the coefficient gradient.
#'
#' # The reordering
#'
#' \pkg{distributions7} keys its Hessian components lexicographically and this
#' package keys them diagonals first, so the second derivatives are looked up
#' under both spellings of each pair, `a_b` and `b_a`, and returned under this
#' package's. The key is composed from the pair rather than parsed out of a
#' name, a parameter name being free to contain an underscore.
#'
#' @param pen A [DistribPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of the parent's free parameters.
#' @param scale Read by the generic, which applies the chain rule onto the
#'   unconstrained scale after dispatch, using the parent's own links. These
#'   methods always return the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` a list of one number per hyperparameter,
#'   named by `pen@params`.
#'   `penalty_hess_theta()` a list of one number per unordered pair, keyed
#'   diagonals first.
#'   `penalty_cross()` a list of one numeric vector of length `pen@n_coef` per
#'   hyperparameter.
#'
#' @examples
#' b <- c(1, 0, -1)
#' pen <- lasso_penalty(n_coef = 3)
#'
#' # The lasso's value is lambda * L1 - q log(lambda/2), so its derivative in
#' # lambda is L1 - q/lambda and the second is q/lambda^2.
#' penalty_grad_theta(pen, b, list(lambda = 0.5))
#' sum(abs(b)) - 3 / 0.5
#' penalty_hess_theta(pen, b, list(lambda = 0.5))
#' 3 / 0.5^2
#'
#' # The mixed block is the sign vector, the value being lambda times the L1
#' # norm in every coefficient direction.
#' penalty_cross(pen, b, list(lambda = 0.5))
#'
#' # An elastic net has two hyperparameters, so three Hessian keys.
#' names(penalty_hess_theta(elasticnet_penalty(n_coef = 3), b,
#'                          list(lambda = 1, alpha = 0.4)))
#'
#' @seealso [penalty_value.DistribPenalty()] for the quantity differentiated,
#'   [penalty_grad_theta()] for the generic and its two scales,
#'   [distributions7::distrib_cross_y()] for the mixed block's source.
#' @keywords internal
S7::method(penalty_grad_theta, DistribPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    g <- distributions7::distrib_gradient(pen@parent, dp_arg(pen, t), theta)
    lapply(g, function(v) -sum(v))
  }

#' @rdname penalty_grad_theta.DistribPenalty
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

#' @rdname penalty_grad_theta.DistribPenalty
#' @name penalty_cross.DistribPenalty
#' @keywords internal
S7::method(penalty_cross, DistribPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    a <- dp_arg(pen, t)
    cy <- distributions7::distrib_cross_y(pen@parent, a, theta)
    lapply(cy, function(v) -map_back(pen, dp_flat(pen, v + 0 * a)))
  }

#' @title Smoothness and Properness of a Separable Penalty
#' @name penalty_kinks.DistribPenalty
#'
#' @description
#' `penalty_kinks()` returns the kinks recorded on the object at construction,
#' and `is_proper()` returns `TRUE` unless the map has fewer rows than columns.
#'
#' @details
#' The kinks do not move with `theta`, unlike SCAD's and MCP's: they are
#' properties of where the parent's location is held, settled by
#' [distrib_kinks()] when the penalty was built. The lasso and the elastic net
#' report `0`; every other shipped parent reports `numeric(0)`.
#'
#' Properness here is a test on the map's shape. The parent is a
#' density by construction, so \eqn{\exp(-\rho)} integrates over the values the
#' penalty reads; what can fail is that a map with fewer rows than columns
#' leaves directions of \eqn{\beta} the penalty never sees, and along those the
#' integral diverges. With no map, or with a map at least as tall as it is
#' wide, the answer is `TRUE`.
#'
#' @param pen A [DistribPenalty()] object.
#' @param theta A named list of the parent's free parameters. Read by
#'   `penalty_kinks()`, whose answer does not depend on it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` a numeric vector, `pen@kinks` as it stands.
#'   `is_proper()` a single logical.
#'
#' @examples
#' # The lasso has a kink at zero and does not move with lambda.
#' penalty_kinks(lasso_penalty(n_coef = 3), list(lambda = 1))
#' penalty_kinks(lasso_penalty(n_coef = 3), list(lambda = 20))
#'
#' # A Gaussian or Student t prior has none.
#' penalty_kinks(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#'
#' # A map that is shorter than it is wide leaves directions unpenalized.
#' c(short = is_proper(lasso_penalty(map = diff(diag(3), differences = 2))),
#'   none  = is_proper(lasso_penalty(n_coef = 3)))
#'
#' @seealso [penalty_kinks()] and [is_proper()] for the generics,
#'   [distrib_kinks()] for where the kinks came from,
#'   [penalty_prox()] for the operator a kink calls for.
#' @keywords internal
S7::method(penalty_kinks, DistribPenalty) <- function(pen, theta, ...) {
  pen@kinks
}

#' @rdname penalty_kinks.DistribPenalty
#' @name is_proper.DistribPenalty
#' @keywords internal
S7::method(is_proper, DistribPenalty) <- function(pen, ...) {
  is.null(pen@map) || nrow(pen@map) >= ncol(pen@map)
}

#' @title What a Separable Penalty's Hyperparameters Are About
#' @name penalty_readable.DistribPenalty
#'
#' @description
#' Returns `NULL` for a univariate parent, whose hyperparameters are already
#' the quantities a reader reads, and
#' [distributions7::mv_derived()] for a multivariate one, whose
#' hyperparameters are the free values of a matrix parameter.
#'
#' @details
#' A univariate parent's hyperparameters are a scale, a rate, a shape, each
#' read on its own scale already, so there is nothing to derive and `NULL` says
#' so. A multivariate parent's are log-Cholesky coordinates or the like, which
#' nobody interprets; what the prior is about is the standard deviations and
#' the correlations of the effects within a block, and the parent declares
#' them.
#'
#' @param pen A [DistribPenalty()] object.
#' @param theta A named list of the parent's free parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `NULL` when `pen@block` is `1L`. Otherwise the list
#'   [distributions7::mv_derived()] returns: `value`, `jacobian`, `transform`
#'   and `block`.
#'
#' @examples
#' # A univariate parent: the rate is the quantity, so nothing is derived.
#' penalty_readable(lasso_penalty(n_coef = 3), list(lambda = 1))
#'
#' # A bivariate Gaussian prior: two standard deviations and a correlation.
#' mv <- distrib_penalty(
#'   distributions7::fixed(distributions7::mvgaussian1_distrib(2),
#'                         mu1 = 0, mu2 = 0), n_coef = 6)
#' r <- penalty_readable(mv, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
#'                                sigma_L2.1 = 0.5))
#' r$value
#' r$transform
#'
#' @seealso [penalty_readable()] for the generic,
#'   [distributions7::mv_derived()] for the declaration this passes through.
#' @keywords internal
S7::method(penalty_readable, DistribPenalty) <- function(pen, theta, ...) {
  # only a blockwise parent has coordinates that are not the quantities: a
  # univariate one's hyperparameters are a scale, a rate, a shape, each read
  # on its own scale already
  if (pen@block == 1L) return(NULL)
  distributions7::mv_derived(pen@parent, theta)
}
