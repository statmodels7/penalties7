#' @include generics.R
NULL

#' @title S7 Class for the Separable Penalty
#'
#' @description
#' `DistribPenalty` is the class of the penalties that [distrib_penalty()]
#' builds. Such a penalty is formed by applying a \pkg{distributions7}
#' log-density to the successive blocks of
#' \eqn{D\beta}. Beyond the seven properties that every penalty carries, the
#' class stores the parent distribution, the points where the parent's
#' log-density has a kink, the block width (one for a univariate parent and
#' the dimension for a multivariate one) and the log-Jacobian of the map.
#'
#' @details
#' The hyperparameters are the parent's own free parameters. The properties
#' `params`, `params_bounds`, `link_params` and `params_smooth` are copied
#' from the distribution object, so the hyperparameters of the penalty, their
#' bounds, their links and their smoothness are those of the parent.
#'
#' @section Methods:
#' The list below gives the methods registered on this class, by the page that
#' documents each.
#'
#' * The value and its coefficient derivatives:
#'   [penalty_value()][penalty_value.DistribPenalty],
#'   [penalty_gradient()][penalty_gradient.DistribPenalty] and
#'   [penalty_hessian()][penalty_hessian.DistribPenalty].
#' * The hyperparameter derivatives:
#'   [penalty_grad_theta()][penalty_grad_theta.DistribPenalty],
#'   [penalty_hess_theta()][penalty_hess_theta.DistribPenalty] and
#'   [penalty_cross()][penalty_cross.DistribPenalty].
#' * The predicates and the structure of the prior:
#'   [penalty_kinks()][penalty_kinks.DistribPenalty],
#'   [is_proper()][is_proper.DistribPenalty],
#'   [penalty_null_basis()][penalty_null_basis.DistribPenalty] and
#'   [penalty_readable()][penalty_readable.DistribPenalty].
#' * The higher derivatives of a marginal criterion:
#'   [penalty_dhessian()][penalty_dhessian.DistribPenalty],
#'   [penalty_d2hessian()][penalty_d2hessian.DistribPenalty],
#'   [penalty_dcross()][penalty_dcross.DistribPenalty] and
#'   [beta_quadratic()][beta_quadratic.DistribPenalty].
#' * The movement of the Hessian with the coefficients:
#'   [penalty_dhessian_beta()][penalty_dhessian_beta.DistribPenalty],
#'   [penalty_d2hessian_beta()][penalty_d2hessian_beta.DistribPenalty] and
#'   [penalty_dhessian_beta_theta()][penalty_dhessian_beta_theta.DistribPenalty].
#' * The proximal operator and the draws:
#'   [penalty_prox()][penalty_prox.DistribPenalty],
#'   [penalty_prox_spec()][penalty_prox_spec.DistribPenalty] and
#'   [penalty_draw()][penalty_draw.DistribPenalty].
#'
#' The class registers no [is_quadratic()][is_quadratic.penalty],
#' [penalty_matrix()][penalty_matrix.penalty],
#' [penalty_rank()][penalty_rank.penalty] or
#' [penalty_logpdet()][penalty_logpdet.penalty] method, so the methods of the
#' base class [penalty()] are used for those generics.
#'
#' @inheritParams penalty
#' @param parent The \pkg{distributions7} object that supplies the
#'   log-density, usually a [distributions7::fixed()] wrapper holding the
#'   location at zero.
#' @param kinks A numeric vector of the points where the parent's log-density
#'   is not differentiable in its argument, possibly empty. [distrib_penalty()]
#'   derives it with [distrib_kinks()] unless its argument `kinks` is supplied.
#' @param block The block width: `1L` for a univariate parent and the parent's
#'   `n_dim` for a multivariate one.
#' @param log_jacobian \eqn{\log\lvert\det D\rvert}, computed by
#'   [map_log_jacobian()][dp_diag_entries] when the penalty is built: `0`
#'   without a map, a finite number for a square nonsingular map, and `NA` for
#'   a map that is not square or is singular. The default is `NA`.
#'
#' @return An S7 object of class `DistribPenalty`, inheriting from [penalty()],
#'   with the seven inherited properties and the four above.
#'
#' @seealso [distrib_penalty()] for the constructor,
#'   [ridge_penalty()] and its siblings for the named instances,
#'   [penalty_value.DistribPenalty()] for the value of the penalty.
#'
#' @examples
#' pen <- lasso_penalty(n_coef = 2)
#' S7::S7_inherits(pen, DistribPenalty)
#'
#' # The parent is a fixed Laplace, the kink is at zero, the blocks are single
#' # coordinates, and the hyperparameter is the rate of the parent.
#' pen@parent
#' pen@kinks
#' pen@block
#' pen@params
#'
#' # Without a map the log-Jacobian is zero.
#' pen@log_jacobian
#'
#' @keywords internal
#' @export
DistribPenalty <- S7::new_class(
  name = "DistribPenalty",
  parent = penalty,
  properties = list(
    parent = S7::class_any,
    kinks = S7::class_numeric,
    block = S7::class_integer,
    log_jacobian = S7::new_property(S7::class_numeric, default = NA_real_)
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
#' # Hyperparameters and derivatives
#'
#' The properties `params`, `params_bounds`, `link_params` and
#' `params_smooth` are copied from the distribution, and every derivative is
#' the distribution's own, carried through the map:
#'
#' \deqn{\frac{\partial\rho}{\partial\beta} = -D'\ell^{(y)}, \qquad
#'   \frac{\partial^2\rho}{\partial\beta^2}
#'     = -D'\,\mathrm{diag}(\ell^{(yy)})\,D, \qquad
#'   \frac{\partial^2\rho}{\partial\beta\,\partial\theta_k}
#'     = -D'\ell^{(y\theta_k)},}
#'
#' where the middle matrix is block diagonal instead of diagonal when the
#' parent is multivariate. The hyperparameter derivatives are the parent's
#' score and Hessian summed over the blocks, and the mixed derivative comes
#' from [distributions7::distrib_cross_y()].
#'
#' # Normalizing constant and Jacobian
#'
#' The normalizing constant of the density is kept, so the value is the
#' negative log prior density. This makes a free scale or a free \eqn{\nu}
#' estimable by a marginal criterion. If the constant were dropped, the value
#' of a Gaussian prior would decrease toward zero as its scale grows, and a
#' criterion would send the scale to infinity.
#'
#' Under a square nonsingular map the value also carries the Jacobian of
#' \eqn{\beta \mapsto D\beta}: it is
#' \eqn{-\sum_i \log f(b_i;\theta) - \log\lvert\det D\rvert}, the negative
#' log-density of \eqn{\beta} instead of that of \eqn{D\beta}. Under a
#' diagonal map, which is the map that standardization builds, the term is
#' \eqn{\sum_j \log\lvert d_j\rvert}. The term is constant in \eqn{\beta} and
#' in \eqn{\theta}, so it changes no derivative, and it makes a marginal
#' criterion invariant to the units of a covariate. The determinant is
#' computed once, when the penalty is built, and stored in the `log_jacobian`
#' property; for a diagonal map it is recomputed from the entries each time
#' the value is computed.
#' Under a map that is not square, or is singular,
#' \eqn{\exp(-\rho)} does not integrate to one over \eqn{\beta},
#' [is_proper()] returns `FALSE`, and the value is
#' \eqn{-\sum_i \log f(b_i;\theta)} with no further constant.
#'
#' A square base matrix or \pkg{Matrix} object whose off-diagonal entries are
#' all exactly zero is stored as a [Matrix::Diagonal()] map. The proximal
#' operator, the draws and the Jacobian recognize a diagonal map by its class,
#' so this conversion makes their result depend on the entries of the map and
#' not on the class in which it was supplied.
#'
#' # A multivariate parent
#'
#' The caller centers the parent, typically through
#' [distributions7::fixed()] at a zero mean, and the matrix parameter carries
#' the dependence within a block. The number of values at which the penalty
#' is read must divide into whole blocks; otherwise the constructor signals an
#' error that names the block width and the count.
#'
#' A blockwise penalty has no proximal operator, because that operator
#' acts one coordinate at a time and the coordinates of a block do not
#' separate. [has_prox()] returns `FALSE` and [penalty_prox()] signals an
#' error.
#'
#' @param d A continuous \pkg{distributions7} object, univariate or
#'   multivariate; typically a [distributions7::fixed()] wrapper holding the
#'   location at zero. A univariate object must inherit from
#'   [distributions7::continuous_distrib]; otherwise the constructor signals
#'   an error that names the distribution. A multivariate object is rejected
#'   in the same way when it is, or wraps, the \pkg{distributions7}
#'   multinomial; a multivariate family written outside \pkg{distributions7}
#'   is not checked. An object whose `dimension` is neither `"univariate"`
#'   nor `"multivariate"` is rejected.
#' @param map The matrix \eqn{D}, with one column per coefficient, or `NULL`
#'   (the default) for the identity. A \pkg{Matrix} object keeps its storage,
#'   so a diagonal or sparse map is not converted to a dense matrix, and a
#'   square map whose off-diagonal entries are all zero is stored as a
#'   [Matrix::Diagonal()].
#' @param n_coef The number of coefficients. It is required when `map` is
#'   `NULL`; otherwise it is ignored, because the count is `ncol(map)`.
#' @param kinks The points where the parent's log-density is not
#'   differentiable in its argument. With `NULL`, the default, they are
#'   derived by [distrib_kinks()] for a univariate parent and set to
#'   `numeric(0)` for a multivariate one, because a kink is a point of a scalar
#'   argument. A numeric vector is used as given, and `numeric(0)` sets the
#'   penalty to have no kinks.
#'
#' @return A [DistribPenalty()] object whose hyperparameters, bounds, links
#'   and `params_smooth` are those of the parent distribution.
#'
#' @examples
#' # A Gaussian prior at zero is the separable form of the ridge, written by
#' # its standard deviation instead of a precision.
#' d <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
#' pen <- distrib_penalty(d, n_coef = 3)
#' pen@params
#' penalty_value(pen, c(1, 0, -1), list(sigma = 2))
#'
#' # The same number comes from stats::dnorm().
#' -sum(stats::dnorm(c(1, 0, -1), sd = 2, log = TRUE))
#'
#' # A correlated prior is built over two coefficients per block and three
#' # blocks. The hyperparameters are the matrix parameter's free values.
#' mv <- distributions7::fixed(distributions7::mvgaussian1_distrib(2),
#'                             mu1 = 0, mu2 = 0)
#' pen2 <- distrib_penalty(mv, n_coef = 6)
#' pen2@block
#' pen2@params
#' penalty_value(pen2, c(1, 0, -1, 0.5, 0.2, -0.3),
#'               list(sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.4))
#'
#' # The blockwise penalty has no proximal operator.
#' has_prox(pen2)
#'
#' # A count that does not divide into blocks is rejected.
#' try(distrib_penalty(mv, n_coef = 5))
#'
#' # A square nonsingular map adds -log|det D| to the value, and a diagonal
#' # base matrix is stored as a Matrix::Diagonal.
#' D <- matrix(c(2, 1, 0, 3), 2)
#' pen3 <- distrib_penalty(d, map = D)
#' penalty_value(pen3, c(0.4, -0.7), list(sigma = 1))
#' -sum(stats::dnorm(D %*% c(0.4, -0.7), log = TRUE)) - log(abs(det(D)))
#' methods::is(distrib_penalty(d, map = diag(c(2, 3)))@map, "diagonalMatrix")
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
    stop(sprintf(paste0(
      "'d' has dimension '%s'; a separable penalty needs a univariate or a\n",
      "  multivariate distribution."), d@dimension), call. = FALSE)
  }
  # the multivariate branch of distributions7 has no continuous/discrete
  # split, so the one discrete family it ships is told by its class, at any
  # depth of wrapping
  discrete <- if (mv) dp_is_multinomial(d) else
    !S7::S7_inherits(d, distributions7::continuous_distrib)
  if (discrete) {
    stop(sprintf(paste0(
      "'d' is '%s', which is not a continuous distribution; a separable\n",
      "  penalty needs a continuous density."), d@distrib_name), call. = FALSE)
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
    # a map whose off-diagonal entries are all zero is a rescaling whatever
    # its class, and is stored as one: the proximal operator, the draws and
    # the Jacobian read a diagonal map by its class
    if (!methods::is(map, "diagonalMatrix") && nrow(map) == ncol(map) &&
        isTRUE(Matrix::isDiagonal(map))) {
      map <- Matrix::Diagonal(x = as.numeric(Matrix::diag(map)))
    }
    q <- ncol(map)
  }
  rows <- if (is.null(map)) q else nrow(map)
  if (rows %% blk != 0L) {
    stop(sprintf(paste0(
      "'%s' is %d-variate, so the %d values at which it is read must divide\n",
      "  into whole blocks; %d does not."),
      d@distrib_name, blk, rows, rows), call. = FALSE)
  }
  DistribPenalty(
    penalty_name = sprintf("separable [%s]", d@distrib_name),
    map = map,
    n_coef = q,
    params = d@params,
    params_bounds = d@params_bounds,
    link_params = d@link_params,
    params_smooth = dp_params_smooth(d),
    parent = d,
    kinks = kinks,
    block = blk,
    log_jacobian = map_log_jacobian(map)
  )
}

#' Helpers of the Separable Constructor
#'
#' @description
#' `dp_params_smooth()` returns the parent's `params_smooth` aligned to its
#' `params`, with `TRUE` for a parameter that `params_smooth` does not list.
#' `dp_is_multinomial()` returns `TRUE` when the object, or any object that
#' it wraps through `parent_distrib`, is a \pkg{distributions7} multinomial.
#'
#' @param d A \pkg{distributions7} object.
#'
#' @return `dp_params_smooth()` returns a named logical vector, one entry per
#'   element of `d@params`. `dp_is_multinomial()` returns a single logical.
#'
#' @seealso [distrib_penalty()]
#'
#' @keywords internal
dp_params_smooth <- function(d) {
  out <- stats::setNames(rep(TRUE, length(d@params)), d@params)
  src <- d@params_smooth
  hit <- intersect(names(src), d@params)
  out[hit] <- as.logical(src[hit])
  out
}

#' @rdname dp_params_smooth
#' @keywords internal
dp_is_multinomial <- function(d) {
  repeat {
    if ("distributions7::MultinomialDistrib" %in% class(d)) return(TRUE)
    if (!"parent_distrib" %in% S7::prop_names(d)) return(FALSE)
    d <- d@parent_distrib
  }
}


#' The Diagonal and the Log-Jacobian of a Map
#'
#' @description
#' `dp_diag_entries()` returns the entries of a \pkg{Matrix} diagonal map,
#' including the zeros, and `NULL` for any other map. `map_log_jacobian()`
#' returns \eqn{\log\lvert\det D\rvert}: `0` for a `NULL` map,
#' \eqn{\sum_j\log\lvert d_j\rvert} for a diagonal one, the log-modulus of the
#' determinant for a square nonsingular one, and `NA` for a map that is not
#' square or is singular. `dp_log_jacobian()` returns the same quantity for a
#' built penalty.
#'
#' @details
#' A diagonal map is treated as singular when an entry is zero or missing. Any
#' other square map is converted once to a dense matrix; it is treated as
#' singular when an entry is not finite or when the rank that [qr()] reports at
#' its default tolerance is below the number of columns, and otherwise the
#' log-modulus comes from [determinant()].
#'
#' `dp_log_jacobian()` recomputes the quantity from the entries of a diagonal
#' map, at the cost of \eqn{q} logarithms, so that the result stays correct
#' when a caller replaces the map with a subset of its diagonal through
#' [S7::set_props()] and the stored property is out of date. For any other
#' map it reads the stored `log_jacobian`.
#'
#' @param map The map \eqn{D}: `NULL`, a base matrix or a \pkg{Matrix}.
#' @param pen A [DistribPenalty()] object.
#'
#' @return `dp_diag_entries()` returns a numeric vector of length `nrow(map)`,
#'   a unit diagonal (`diag = "U"`) giving ones, or `NULL`. `map_log_jacobian()`
#'   and `dp_log_jacobian()` return a single number, possibly `NA`.
#'
#' @seealso [distrib_penalty()], [penalty_value.DistribPenalty()]
#'
#' @keywords internal
dp_diag_entries <- function(map) {
  if (is.null(map) || !isS4(map) || !methods::is(map, "diagonalMatrix")) {
    return(NULL)
  }
  if (identical(map@diag, "U")) rep(1, nrow(map)) else as.numeric(map@x)
}

#' @rdname dp_diag_entries
#' @keywords internal
map_log_jacobian <- function(map) {
  if (is.null(map)) return(0)
  d <- dp_diag_entries(map)
  if (!is.null(d)) {
    if (anyNA(d) || any(d == 0)) return(NA_real_)
    return(sum(log(abs(d))))
  }
  if (nrow(map) != ncol(map)) return(NA_real_)
  m <- as.matrix(map)
  if (!all(is.finite(m)) || qr(m)$rank < ncol(m)) return(NA_real_)
  as.numeric(determinant(m, logarithm = TRUE)$modulus)
}

#' @rdname dp_diag_entries
#' @keywords internal
dp_log_jacobian <- function(pen) {
  if (is.null(pen@map)) return(0)
  if (!is.null(dp_diag_entries(pen@map))) return(map_log_jacobian(pen@map))
  pen@log_jacobian
}


#' Where the Parent's Log-Density Has a Kink
#'
#' @description
#' Returns the points at which a separable penalty built on this distribution
#' is not differentiable, derived from the distribution itself, or
#' `numeric(0)` when there are none. Among the families that
#' \pkg{distributions7} ships, only the two Laplace parametrizations and the
#' elastic net have a location in which the log-density is not smooth, so only
#' these give a kink, and only when [distributions7::fixed()] holds the
#' location.
#'
#' @details
#' # How a candidate is found
#'
#' A distribution records in `params_smooth` the parameters in which its
#' log-likelihood is differentiable, and for a location family a location
#' that is not smooth is a kink in the argument at that location. A penalty
#' is the negative log-density read at the coefficient, so a
#' [distributions7::fixed()] wrapper that holds such a parameter at a value
#' puts the kink at that value: `fixed(laplace_distrib(), mu = 0)` is a
#' Laplace prior at zero and has a kink at zero, while
#' `fixed(laplace_distrib(), mu = 2)` has one at two.
#'
#' A free parameter gives no candidate, because its value is set by the
#' hyperparameters at the time of each evaluation. A Laplace with a free
#' location therefore has no reported kink.
#'
#' The wrapper is recognized by the properties that it carries,
#' `parent_distrib` and `fixed_params`, and not by its class, because the
#' classes of the [distributions7::fixed()] wrapper are not exported and a
#' family written outside \pkg{distributions7} may hold its parameters in the
#' same way. The chain of `parent_distrib` is walked from the outermost object
#' down, and every level that carries `fixed_params` contributes candidates,
#' so that `truncated(fixed(laplace2_distrib(), mu = 0), lower = -1,
#' upper = 3)` has its kink at zero.
#'
#' # How a candidate is confirmed
#'
#' A candidate is dropped when it is not finite or when it does not lie
#' strictly inside the support of the outermost object, because a point at the
#' edge of the support is not a kink of the density on its support: a
#' `folded()` Laplace held at zero has no kink at zero. Each remaining
#' candidate is then tested by comparing the one-sided derivatives of the
#' log-density across it at a probe point ([probe_theta()], [has_jump()]),
#' and a candidate at which the derivative does not jump is dropped. The
#' inference from `params_smooth` supplies the points to test, and the test
#' removes the candidates at which the density is smooth, for example a
#' candidate that comes from a non-smooth parameter that is not a location.
#' The test is made on the outermost object, so a wrapper that moves the
#' point drops the candidate.
#'
#' @param d A univariate \pkg{distributions7} object. An object without
#'   `parent_distrib`, or whose chain of parents has no level carrying
#'   `fixed_params`, yields `numeric(0)` without further inspection.
#'
#' @return A numeric vector of the confirmed kink points, without duplicates
#'   and possibly empty.
#'
#' @seealso [distrib_penalty()], which calls this function when its argument
#'   `kinks` is `NULL`, [penalty_kinks()] for the kinks of a built penalty,
#'   [has_jump()] for the test.
#'
#' @examples
#' # A Laplace held at a location has a kink at that location.
#' distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
#'                                     mu = 0))
#' distrib_kinks(distributions7::fixed(distributions7::laplace_distrib(),
#'                                     mu = 2))
#'
#' # A Gaussian is smooth everywhere, so there is no candidate.
#' distrib_kinks(distributions7::fixed(distributions7::gaussian1_distrib(),
#'                                     mu = 0))
#'
#' # A free location gives no candidate.
#' distrib_kinks(distributions7::laplace_distrib())
#'
#' # A fixed() level under another wrapper is found; the origin of a folded
#' # Laplace is the edge of its support and is dropped.
#' lap0 <- distributions7::fixed(distributions7::laplace2_distrib(), mu = 0)
#' distrib_kinks(distributions7::truncated(lap0, lower = -1, upper = 3))
#' distrib_kinks(distributions7::folded(lap0))
#'
#' @export
distrib_kinks <- function(d) {
  # the wrapper is recognized by what it carries rather than by its class:
  # the two fixed classes are not exported, and a family written outside
  # distributions7 may hold its parameters the same way. It may sit under
  # other wrappers, so the chain of parents is walked and every fixed level
  # contributes its candidates
  cand <- numeric(0)
  inner <- d
  while ("parent_distrib" %in% S7::prop_names(inner)) {
    if ("fixed_params" %in% S7::prop_names(inner)) {
      smooth <- inner@parent_distrib@params_smooth
      fixed_at <- inner@fixed_params
      rough <- names(smooth)[!smooth]
      cand <- c(cand, unlist(fixed_at[intersect(rough, names(fixed_at))],
                             use.names = FALSE))
    }
    inner <- inner@parent_distrib
  }
  # a point at the edge of the outer support (a folded parent's origin) is
  # not a kink of the density on its support
  b <- d@bounds
  cand <- unique(cand[is.finite(cand) & cand > b[1L] & cand < b[2L]])
  if (!length(cand)) return(numeric(0))
  # confirmed on the OUTER object: a wrapper that moves the point drops the
  # candidate here
  th <- probe_theta(d)
  as.numeric(cand[vapply(cand, function(k) has_jump(d, th, k), logical(1))])
}


#' A Point at Which to Read the Parent
#'
#' @description
#' Returns one value per free parameter, at which the parent can be evaluated
#' when the hyperparameters are not given. [distrib_kinks()] uses it to
#' differentiate the log-density before any penalty exists.
#'
#' @details
#' The value of a parameter is the midpoint of its bounds where both are
#' finite, one above the lower bound where only that bound is finite, one
#' below the upper bound where only that bound is finite, and `1` where
#' both bounds are infinite. A scale on \eqn{(0, \infty)} is therefore probed at
#' 1 and a mixing weight on \eqn{(0, 1)} at 0.5.
#'
#' The rule assumes that whether a point is a kink does not depend on the
#' values of the other parameters, so that any admissible point serves.
#'
#' @param d A \pkg{distributions7} object.
#'
#' @return A named list of one number per element of `d@params`, in that order.
#'   The value is `list()` for a distribution with no free parameters.
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


#' Test for a Jump in the Slope of the Log-Density
#'
#' @description
#' Compares the derivative of the log-density in its argument on either side
#' of a candidate point at two widths, and returns `TRUE` when the difference
#' has the behavior of a jump. [distrib_kinks()] uses it to confirm a
#' candidate kink inferred from `params_smooth`.
#'
#' @details
#' With \eqn{s = \max(1, \lvert a\rvert)}, where \eqn{a} is `at`, the
#' difference \eqn{g(h) = \lvert\ell^{(y)}(a + h) - \ell^{(y)}(a - h)\rvert}
#' is computed at the narrow width \eqn{h = \varepsilon s} and at the wide
#' width \eqn{h = 10\varepsilon s}. At a kink the difference is of the order
#' of the jump at both widths (twice the rate for a Laplace at its location)
#' and does not grow with the width. At a smooth point it is about
#' \eqn{2h\lvert\ell^{(yy)}\rvert}, which can exceed any fixed threshold when
#' the curvature is large, and it grows tenfold from the narrow width to the
#' wide one. A kink is therefore declared when \eqn{g(\varepsilon s)} exceeds
#' \eqn{10^{-6}s}, a threshold that is absolute near the origin and relative
#' away from it, and is more than half of \eqn{g(10\varepsilon s)}. The
#' widths are scaled up from \eqn{\varepsilon} and not down, so that a slope
#' that the parent obtains by its own numerical differentiation still
#' resolves the jump at the narrow width.
#'
#' When [distributions7::distrib_grad_y()] returns a value that is not finite
#' or signals an error at any of the four points, the result is `FALSE`, so
#' the candidate is dropped.
#'
#' @param d A \pkg{distributions7} object.
#' @param theta A named list of parameter values at which the log-density is
#'   evaluated, as [probe_theta()] returns.
#' @param at The candidate point, a single number.
#' @param eps The relative narrow width, `1e-5` by default; the wide width is
#'   ten times as large, and both are multiplied by \eqn{\max(1,
#'   \lvert a\rvert)}. It should be large enough that the one-sided
#'   derivatives are not dominated by rounding, and small enough that a
#'   second feature of the density does not fall inside the window.
#'
#' @return A single logical, `TRUE` when the slope jumps.
#'
#' @seealso [distrib_kinks()], [probe_theta()]
#'
#' @keywords internal
has_jump <- function(d, theta, at, eps = 1e-5) {
  # at a kink the jump does not grow with the width; at a smooth point it is
  # about 2 h |l''| and grows tenfold with h. The widths are widened rather
  # than narrowed so that a slope differenced by the parent's own fallback
  # still resolves the jump at the narrower one
  s <- max(1, abs(at))
  g <- tryCatch(vapply(c(1, 10) * eps * s, function(h) {
    up <- distributions7::distrib_grad_y(d, at + h, theta)
    dn <- distributions7::distrib_grad_y(d, at - h, theta)
    abs(as.numeric(up) - as.numeric(dn))
  }, numeric(1)), error = function(e) c(NA_real_, NA_real_))
  isTRUE(all(is.finite(g)) && g[1L] > 1e-6 * s && g[1L] > g[2L] / 2)
}


#' @title The Named Penalties
#'
#' @description
#' These functions construct the four penalties that a modeling layer names
#' directly: `ridge_penalty()`, `lasso_penalty()`, `elasticnet_penalty()` and
#' `heavy_penalty()`. Each is a prior centered at zero. In the first three the
#' hyperparameter `lambda` is a precision or a rate, and a larger value
#' shrinks harder. In `heavy_penalty()` the hyperparameter `sigma` is a scale,
#' and a larger value shrinks less.
#'
#' @details
#' # Ridge
#'
#' `ridge_penalty()` is the Gaussian prior at zero written by its
#' precision, which is [quadratic_penalty()] at the identity matrix: its
#' value agrees with that of the separable Gaussian prior at
#' \eqn{\sigma = \lambda^{-1/2}} up to rounding. It is therefore built on the
#' quadratic branch, and its hyperparameter is the `lambda` of that branch. It
#' is the only one of the four that is not a `DistribPenalty`, and the only
#' one for which [is_quadratic()] returns `TRUE`.
#'
#' The separable form is
#' `distrib_penalty(fixed(gaussian1_distrib(), mu = 0), n_coef = q)`, whose
#' hyperparameter is the standard deviation \eqn{\sigma} instead. It is the
#' form to use when a standard deviation is the quantity to be reported.
#'
#' # Lasso
#'
#' The lasso is the Laplace in location and rate (`laplace2`) held at zero, so
#' the free hyperparameter is the rate \eqn{\lambda} and
#'
#' \deqn{\rho(\beta;\lambda) = \lambda\lVert D\beta\rVert_1
#'   - m\log(\lambda/2),}
#'
#' where \eqn{m} is the number of values at which the penalty is read (the
#' number of rows of \eqn{D}, or the number of coefficients without a map).
#' Under a square nonsingular map the value also carries
#' \eqn{-\log\lvert\det D\rvert}, as described in [distrib_penalty()]. The
#' kink at zero is declared, so [penalty_kinks()] reports it and
#' [penalty_prox()] gives the soft threshold.
#'
#' # Elastic net
#'
#' The elastic net is the product of the Laplace and the Gaussian at zero,
#' normalized ([distributions7::enet_distrib()]), so the hyperparameters are
#' the overall rate \eqn{\lambda} and the mixing weight \eqn{\alpha} on \eqn{(0,
#' 1)}:
#'
#' \deqn{\rho(\beta;\lambda,\alpha) = \lambda\left\{
#'   \alpha\lVert D\beta\rVert_1
#'   + (1-\alpha)\lVert D\beta\rVert_2^2/2\right\} + m\,c(\lambda, \alpha),}
#'
#' where \eqn{m} is as above and
#' \eqn{c(\lambda,\alpha) = \log\int \exp[-\lambda\{\alpha\lvert x\rvert +
#' (1-\alpha)x^2/2\}]\,dx} is the constant of one coordinate. The same
#' Jacobian term is added under a square nonsingular map. The constant
#' depends on both hyperparameters, so a marginal criterion can estimate them.
#' Without the constant, both hyperparameters would be unidentified.
#'
#' # Heavy-tailed prior
#'
#' The heavy-tailed prior is the Student t at zero, with a free scale
#' \eqn{\sigma} and free degrees of freedom \eqn{\nu}. The normalizing constant
#' is kept, which makes \eqn{\nu} estimable by a marginal criterion. With a
#' small \eqn{\nu} the prior shrinks small coefficients as a Gaussian prior does
#' and leaves large ones nearly unchanged, which a Gaussian prior does not do at
#' any scale. Like the ridge, it has no kink, so it does not produce exact
#' zeros.
#'
#' @param map The matrix \eqn{D}, with one column per coefficient, or `NULL`
#'   (the default) for the identity. If a map is given, the number of
#'   coefficients is `ncol(map)` and `n_coef` is ignored.
#' @param n_coef The number of coefficients when `map` is `NULL`. It is a
#'   single whole number, and the default is `1L`.
#'
#' @return `ridge_penalty()` returns a [QuadraticPenalty()] with the
#'   hyperparameter `lambda`, a precision, on \eqn{(0, \infty)}.
#'   `lasso_penalty()` returns a [DistribPenalty()] with `lambda`, a rate, on
#'   \eqn{(0, \infty)} and a kink at zero. `elasticnet_penalty()` returns a
#'   [DistribPenalty()] with `lambda` on \eqn{(0, \infty)} and `alpha` on
#'   \eqn{(0, 1)}, and a kink at zero. `heavy_penalty()` returns a
#'   [DistribPenalty()] with `sigma` and `nu` both on \eqn{(0, \infty)}, and no
#'   kink.
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
#' # The ridge is on the quadratic branch, and its lambda is a precision.
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
#' # Of these four, only the lasso and the elastic net have a kink, so only
#' # they can set a coefficient exactly to zero.
#' penalty_kinks(l, list(lambda = 1))
#' penalty_kinks(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#'
#' # The heavy-tailed prior shrinks the large coefficient far less than a
#' # unit ridge, and the small ones somewhat more.
#' penalty_prox(heavy_penalty(n_coef = 3), c(4, 0.3, -1), 1,
#'              list(sigma = 1, nu = 3))
#' penalty_prox(r, c(4, 0.3, -1), 1, list(lambda = 1))
#'
#' @seealso [distrib_penalty()] for the construction of the lasso, the
#'   elastic net and the heavy-tailed prior,
#'   [quadratic_penalty()] for the construction of the ridge,
#'   [scad_penalty()] and [mcp_penalty()] for the non-convex alternatives,
#'   [penalty_prox()] for their proximal operators.
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
#' `dp_arg()` reshapes \eqn{D\beta} into the argument of the parent: the
#' vector itself for a univariate parent, and a matrix of one row per block for
#' a multivariate one. `dp_flat()` undoes the reshaping, carrying the result
#' of the parent back to a vector in coefficient order.
#'
#' @details
#' The blocks are the successive stretches of \eqn{D\beta}, so the
#' reshaping fills by row: block \eqn{i} occupies positions
#' \eqn{(i-1)p+1, \dots, ip}. This is the order in which a grouped design
#' assembles its coefficients and the order that \eqn{I_m \otimes \Sigma}
#' assumes, so the two agree.
#'
#' Both functions return their argument unchanged when the block width is
#' one, which is the case of every univariate parent.
#'
#' @param pen A [DistribPenalty()] object.
#' @param t The mapped coefficient vector \eqn{D\beta}, of length a multiple of
#'   `pen@block`. This argument belongs to `dp_arg()` only.
#' @param g The result of the parent: a vector for a univariate parent, a
#'   matrix of one row per block for a multivariate one. This argument belongs
#'   to `dp_flat()` only.
#'
#' @return `dp_arg()` returns the vector unchanged when `pen@block` is `1L`, and
#'   otherwise a matrix with `pen@block` columns filled by row.
#'   `dp_flat()` returns a numeric vector, in coefficient order.
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
#' block-diagonal middle matrix of the Hessian of a multivariate separable
#' penalty: one \eqn{p \times p} block per block of coefficients, and zeros
#' between them.
#'
#' @details
#' The zeros between the blocks express the separability. Coordinates in
#' different blocks are independent under the prior, so the Hessian has no
#' entry linking them, however dependent the coordinates within a block are.
#'
#' Two input shapes are accepted because \pkg{distributions7} returns two. A
#' family whose response Hessian does not depend on the observation, such as
#' the multivariate Gaussian, returns one \eqn{p \times p} matrix, which is
#' placed in every block; any other family returns a
#' \eqn{p \times p \times n} array, and each slice goes to its own block. The
#' two shapes are distinguished with `is.matrix()`.
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
#' parent summed over the blocks \eqn{b_i} of \eqn{D\beta}. Because the parent
#' supplies its own normalizing constant, the value is the negative log prior
#' density of \eqn{D\beta}. Under a square nonsingular map, diagonal or not,
#' \eqn{\log\lvert\det D\rvert} is subtracted, which makes the value the
#' negative log-density of \eqn{\beta}. Under a map that is not square, or is
#' singular, the function does not subtract a constant.
#'
#' @details
#' The map is applied, the result is reshaped into the argument of the parent
#' by [dp_arg()], and one call to [distributions7::distrib_pdf()] with
#' `log = TRUE` gives the contribution of every block at once. A univariate
#' parent receives the vector, and a \eqn{p}-variate one a matrix of one row
#' per block.
#'
#' The log-determinant is read by [dp_log_jacobian()][dp_diag_entries]: it is
#' computed from the entries of a diagonal map, and for any other map it is
#' the `log_jacobian` stored when the penalty was built.
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
#' # A Gaussian prior at zero: the value is its negative log-density.
#' g0 <- distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0)
#' pen <- distrib_penalty(g0, n_coef = 3)
#' b <- c(1, 0, -1)
#' penalty_value(pen, b, list(sigma = 2))
#' -sum(stats::dnorm(b, sd = 2, log = TRUE))
#'
#' # The value of the lasso is lambda times the L1 norm minus q log(lambda/2).
#' penalty_value(lasso_penalty(n_coef = 3), b, list(lambda = 0.5))
#' 0.5 * sum(abs(b)) - 3 * log(0.5 / 2)
#'
#' # A diagonal map subtracts the sum of the log-moduli of its entries, and a
#' # map with fewer rows than columns subtracts nothing.
#' dg <- distrib_penalty(g0, map = Matrix::Diagonal(x = c(2, 0.5, 4)))
#' penalty_value(dg, b, list(sigma = 2))
#' -sum(stats::dnorm(c(2, 0, -4), sd = 2, log = TRUE)) - log(2 * 0.5 * 4)
#' sh <- distrib_penalty(g0, map = diff(diag(3)))
#' penalty_value(sh, b, list(sigma = 2))
#' -sum(stats::dnorm(diff(b), sd = 2, log = TRUE))
#'
#' @seealso [distrib_penalty()] for the construction,
#'   [penalty_gradient.DistribPenalty()] for the coefficient derivatives,
#'   [ridge_penalty()] and its siblings for the named instances.
#' @keywords internal
S7::method(penalty_value, DistribPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  v <- -sum(distributions7::distrib_pdf(pen@parent, dp_arg(pen, t), theta,
                                        log = TRUE))
  if (is.null(pen@map)) return(v)
  lj <- dp_log_jacobian(pen)
  if (is.na(lj)) v else v - lj
}

#' @title Coefficient Derivatives of a Separable Penalty
#' @name penalty_gradient.DistribPenalty
#'
#' @description
#' These methods return the response derivatives of the parent, negated and
#' carried back through the map. `penalty_gradient()` returns
#' \eqn{-D'\ell^{(y)}} and `penalty_hessian()` returns
#' \eqn{-D'\,\mathrm{diag}(\ell^{(yy)})\,D}, where the middle matrix is block
#' diagonal instead of diagonal when the parent is multivariate.
#'
#' @details
#' The derivatives come from [distributions7::distrib_grad_y()] and
#' [distributions7::distrib_hess_y()], which are in closed form for the
#' continuous families that \pkg{distributions7} ships; for a family that
#' supplies only its density, \pkg{distributions7} differentiates
#' numerically. These methods do not compute any derivative themselves.
#'
#' The `+ 0 * a` in both bodies is a recycling guard: if a family returns a
#' single number where a vector is expected, adding zero times the argument
#' widens it to the length of the argument without changing a value.
#'
#' At a kink the gradient is one element of the subdifferential, chosen by the
#' parent. For the lasso at \eqn{\beta_j = 0} that element is `0`, while the
#' subdifferential is \eqn{[-\lambda, \lambda]}; the exact zero of a fit is
#' produced by [penalty_prox()].
#'
#' @param pen A [DistribPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of the parent's free parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` returns a numeric vector of length `pen@n_coef`;
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
#' # A multivariate parent gives a block-diagonal Hessian, because the
#' # coordinates depend on each other within a block and are independent
#' # between blocks.
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
#' These methods return the parent's own score and information in its
#' parameters, negated and summed over the blocks, plus the mixed block.
#' `penalty_grad_theta()` returns \eqn{-\sum_i \ell^{(\theta_k)}},
#' `penalty_hess_theta()` returns \eqn{-\sum_i \ell^{(\theta_k\theta_l)}}, and
#' `penalty_cross()` returns \eqn{-D'\ell^{(y\theta_k)}}, one coefficient
#' vector per hyperparameter.
#'
#' @details
#' The first two are [distributions7::distrib_gradient()] and
#' [distributions7::distrib_hessian()] evaluated at the blocks and summed. The
#' map does not enter them, because a hyperparameter of the prior does not act
#' through \eqn{D}. The third is [distributions7::distrib_cross_y()], carried
#' back through the map like the coefficient gradient.
#'
#' # Keys of the second derivatives
#'
#' The second derivatives are returned under the keys of [ptheta_pairs()]:
#' the diagonal pairs first, then the off-diagonal pairs, in the order of
#' `pen@params`. The shipped continuous families of \pkg{distributions7} use
#' the same order, but a family may spell a pair the other way, so each pair
#' is looked up under both spellings, `a_b` and `b_a`. The key is composed
#' from the pair and not parsed out of a name, because a parameter name may
#' contain an underscore.
#'
#' @param pen A [DistribPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of the parent's free parameters.
#' @param scale Read by the generic, which applies the chain rule to the
#'   unconstrained scale after dispatch and uses the links copied from the
#'   parent. These methods always return derivatives on the parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` returns a list of one number per
#'   hyperparameter, named by `pen@params`. `penalty_hess_theta()` returns a
#'   list of one number per unordered pair, keyed diagonals first.
#'   `penalty_cross()` returns a list of one numeric vector of length
#'   `pen@n_coef` per hyperparameter.
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
#' # The mixed block is the sign vector, because the value is lambda times the
#' # L1 norm in every coefficient direction.
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
#' and `is_proper()` returns `TRUE` when there is no map or the map is square
#' and nonsingular, and `FALSE` otherwise.
#'
#' @details
#' The kinks do not move with `theta`, unlike those of SCAD and MCP: they are
#' fixed by the value at which the location of the parent is held, and they
#' are set when the penalty is built, by [distrib_kinks()] or by the argument
#' `kinks` of [distrib_penalty()]. Of the named penalties, the lasso and the
#' elastic net report `0` and [heavy_penalty()] reports `numeric(0)`.
#'
#' Properness is a test on the map. The parent is a density, so
#' \eqn{\exp(-\rho)} integrates to one over the values at which the penalty is
#' read, and under a square nonsingular map the value carries
#' \eqn{-\log\lvert\det D\rvert}, so that \eqn{\exp(-\rho)} integrates to one
#' over \eqn{\beta} as well. A map whose rank is below its number of columns
#' leaves directions of \eqn{\beta} on which the penalty does not depend, and
#' along those directions the integral diverges. Under a map with more rows
#' than columns and full column rank the integral does not equal one, and it
#' depends on the hyperparameters through a constant that the value does not
#' carry. The test reads the `log_jacobian` stored at construction, or the
#' entries of a diagonal map, through [dp_log_jacobian()][dp_diag_entries],
#' and returns `TRUE` exactly when that quantity is finite.
#'
#' @param pen A [DistribPenalty()] object.
#' @param theta A named list of the parent's free parameters. It is accepted
#'   by `penalty_kinks()`, whose result does not depend on it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` returns a numeric vector, `pen@kinks` as it stands.
#'   `is_proper()` returns a single logical.
#'
#' @examples
#' # The kink of the lasso is at zero for every lambda.
#' penalty_kinks(lasso_penalty(n_coef = 3), list(lambda = 1))
#' penalty_kinks(lasso_penalty(n_coef = 3), list(lambda = 20))
#'
#' # A Student t prior has none.
#' penalty_kinks(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#'
#' # The penalty is proper without a map and with a square nonsingular map. It
#' # is improper with a map that has more rows than columns, with a singular
#' # one, and with one that has fewer rows than columns.
#' c(none     = is_proper(lasso_penalty(n_coef = 3)),
#'   square   = is_proper(lasso_penalty(map = matrix(c(2, 1, 0, 3), 2))),
#'   tall     = is_proper(lasso_penalty(map = rbind(diag(2), c(1, 1)))),
#'   singular = is_proper(lasso_penalty(map = matrix(1, 2, 2))),
#'   short    = is_proper(lasso_penalty(map = diff(diag(3), differences = 2))))
#'
#' @seealso [penalty_kinks()] and [is_proper()] for the generics,
#'   [distrib_kinks()] for the derivation of the kinks,
#'   [penalty_null_basis.DistribPenalty()] for the directions that an
#'   improper map leaves free,
#'   [penalty_prox()] for the proximal operator at a kink.
#' @keywords internal
S7::method(penalty_kinks, DistribPenalty) <- function(pen, theta, ...) {
  pen@kinks
}

#' @rdname penalty_kinks.DistribPenalty
#' @name is_proper.DistribPenalty
#' @keywords internal
S7::method(is_proper, DistribPenalty) <- function(pen, ...) {
  isTRUE(is.finite(dp_log_jacobian(pen)))
}

#' @title Null Basis of a Separable Penalty
#' @name penalty_null_basis.DistribPenalty
#'
#' @description
#' Returns an orthonormal basis of the null space of the map \eqn{D}, that is,
#' of the directions of \eqn{\beta} on which the penalty does not depend. The
#' basis has no columns when there is no map or when the map has full column
#' rank. Under a diagonal map the basis is made of the unit vectors at the
#' entries of the diagonal that are zero or missing, so it is spanned by
#' coordinates.
#'
#' @details
#' A separable penalty is not a quadratic form, and its log-density is not
#' computed through a pseudo-determinant. A marginal criterion needs from an
#' improper separable prior only the directions that the prior leaves free,
#' and this basis supplies them. The rank of a non-diagonal map is the rank
#' that [qr()] reports at its default tolerance, and the basis is formed by the
#' right singular vectors of the map beyond that rank. For a square map this
#' is the same rank test that [is_proper()] uses. A map with more rows than
#' columns and full column rank has an empty basis, although [is_proper()]
#' returns `FALSE` for it.
#'
#' @param pen A [DistribPenalty()] object.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `pen@n_coef` rows and one column per null
#'   direction.
#'
#' @examples
#' # The first call has no map, and the second has a diagonal map with a zero
#' # at the second coordinate.
#' penalty_null_basis(lasso_penalty(n_coef = 3))
#' penalty_null_basis(lasso_penalty(map = Matrix::Diagonal(x = c(1, 0, 2))))
#'
#' # A map with one row and three columns leaves two directions free.
#' penalty_null_basis(heavy_penalty(map = diff(diag(3), differences = 2)))
#'
#' @seealso [penalty_null_basis()] for the generic, [is_proper()] for the
#'   test of properness.
#' @keywords internal
S7::method(penalty_null_basis, DistribPenalty) <- function(pen, ...) {
  q <- as.integer(pen@n_coef)
  if (is.null(pen@map)) return(matrix(0, q, 0L))
  d <- dp_diag_entries(pen@map)
  if (!is.null(d)) {
    zero <- which(is.na(d) | d == 0)
    return(diag(1, q)[, zero, drop = FALSE])
  }
  m <- as.matrix(pen@map)
  r <- qr(m)$rank
  if (r >= q) return(matrix(0, q, 0L))
  svd(m, nu = 0L, nv = q)$v[, seq.int(r + 1L, q), drop = FALSE]
}

#' @title Interpretable Quantities of a Separable Penalty
#' @name penalty_readable.DistribPenalty
#'
#' @description
#' Returns `NULL` for a univariate parent, whose hyperparameters are already
#' interpretable quantities, and the result of
#' [distributions7::mv_derived()] for a multivariate one, whose
#' hyperparameters are the free values of a matrix parameter.
#'
#' @details
#' The hyperparameters of a univariate parent are quantities such as a scale,
#' a rate or a shape, each on its own scale, so there is nothing to derive and
#' the method returns `NULL`. Those of a multivariate parent are coordinates
#' such as log-Cholesky entries, which are not directly interpretable; the
#' interpretable quantities are the standard deviations and the correlations
#' of the effects within a block, and the parent provides them through
#' [distributions7::mv_derived()].
#'
#' @param pen A [DistribPenalty()] object.
#' @param theta A named list of the parent's free parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `NULL` when `pen@block` is `1L`. Otherwise the method returns the
#'   list that [distributions7::mv_derived()] produces, with the elements
#'   `value`, `jacobian`, `transform` and `block`.
#'
#' @examples
#' # For a univariate parent the rate is already the quantity of interest, so
#' # the result is NULL.
#' penalty_readable(lasso_penalty(n_coef = 3), list(lambda = 1))
#'
#' # A bivariate Gaussian prior has two standard deviations and a correlation.
#' mv <- distrib_penalty(
#'   distributions7::fixed(distributions7::mvgaussian1_distrib(2),
#'                         mu1 = 0, mu2 = 0), n_coef = 6)
#' r <- penalty_readable(mv, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
#'                                sigma_L2.1 = 0.5))
#' r$value
#' r$transform
#'
#' @seealso [penalty_readable()] for the generic,
#'   [distributions7::mv_derived()] for the quantities that a multivariate
#'   family provides.
#' @keywords internal
S7::method(penalty_readable, DistribPenalty) <- function(pen, theta, ...) {
  # only a blockwise parent has coordinates that are not the quantities: a
  # univariate one's hyperparameters are a scale, a rate, a shape, each read
  # on its own scale already
  if (pen@block == 1L) return(NULL)
  distributions7::mv_derived(pen@parent, theta)
}
