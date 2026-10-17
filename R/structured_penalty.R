#' @include quadratic_penalty.R
NULL

#' @title S7 Class for the Structured Quadratic Penalty
#'
#' @description
#' `StructuredPenalty` is the class of the objects that [structured_penalty()]
#' builds: the Gaussian prior whose precision is a \pkg{parameters7} matrix
#' parameter. It adds one property to [penalty()], the structure itself, and
#' every quantity of the branch is computed from the methods of the structure.
#'
#' @details
#' The hyperparameters are the structure's free vector, so `params` is its
#' `free_names`, every entry of `params_bounds` is
#' \eqn{(-\infty, \infty)} and every link is the identity. The constraint that
#' keeps the matrix positive definite is part of the structure, where a scalar
#' link cannot express it. `map` is always `NULL`.
#'
#' @section Methods:
#' This section lists the methods registered on this class, by the page that
#' documents each.
#'
#' * The value and its coefficient derivatives:
#'   [penalty_value()][penalty_value.StructuredPenalty],
#'   [penalty_gradient()][penalty_gradient.StructuredPenalty] and
#'   [penalty_hessian()][penalty_hessian.StructuredPenalty].
#' * The hyperparameter derivatives:
#'   [penalty_grad_theta()][penalty_grad_theta.StructuredPenalty],
#'   [penalty_hess_theta()][penalty_hess_theta.StructuredPenalty] and
#'   [penalty_cross()][penalty_cross.StructuredPenalty].
#' * The predicates: [penalty_kinks()][penalty_kinks.StructuredPenalty],
#'   [is_proper()][is_proper.StructuredPenalty] and
#'   [is_quadratic()][is_quadratic.StructuredPenalty].
#' * The quantities of a marginal criterion:
#'   [penalty_matrix()][penalty_matrix.StructuredPenalty],
#'   [penalty_rank()][penalty_rank.StructuredPenalty],
#'   [penalty_null_basis()][penalty_null_basis.StructuredPenalty] and
#'   [penalty_logpdet()][penalty_logpdet.StructuredPenalty].
#' * The higher derivatives of a marginal criterion:
#'   [penalty_dhessian()][penalty_dhessian.StructuredPenalty],
#'   [penalty_d2hessian()][penalty_d2hessian.StructuredPenalty],
#'   [penalty_dcross()][penalty_dcross.StructuredPenalty] and
#'   [beta_quadratic()][beta_quadratic.StructuredPenalty].
#' * The movement of the Hessian with the coefficients:
#'   [penalty_dhessian_beta()][penalty_dhessian_beta.StructuredPenalty],
#'   [penalty_d2hessian_beta()][penalty_d2hessian_beta.penalty] and
#'   [penalty_dhessian_beta_theta()][penalty_dhessian_beta_theta.penalty].
#'   The last two return zero matrices because the Hessian does not depend on
#'   the coefficients, and are documented on the page of the base-class
#'   methods.
#' * The proximal operator:
#'   [penalty_prox()][penalty_prox.StructuredPenalty], documented on the page
#'   of the quadratic branch.
#' * The interpretable quantities:
#'   [penalty_readable()][penalty_readable.StructuredPenalty], which returns
#'   those that the structure declares through
#'   [parameters7::param_readable()].
#'
#' The class has no [penalty_prox_spec()] or [penalty_draw()] method of its own,
#' so those generics dispatch to the methods of the base class [penalty()].
#'
#' @inheritParams penalty
#' @param structure A \pkg{parameters7} `matrix_parameter`, read as the prior's
#'   precision. It supplies the dimension, the rank, the null basis, the free
#'   names, the matrix and its derivative arrays.
#'
#' @return An S7 object of class `StructuredPenalty`, inheriting from
#'   [penalty()], with the seven inherited properties and `structure`.
#'
#' @seealso [structured_penalty()] for the constructor to use,
#'   [penalty_value.StructuredPenalty()] for what the branch computes,
#'   [quadratic_penalty()] for the branch whose matrix is fixed.
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4))
#' S7::S7_inherits(pen, StructuredPenalty)
#'
#' # The hyperparameters are the structure's free names, unconstrained and
#' # on identity links.
#' pen@params
#' pen@params_bounds
#'
#' @keywords internal
#' @export
StructuredPenalty <- S7::new_class(
  name = "StructuredPenalty",
  parent = penalty,
  properties = list(
    structure = S7::class_any
  )
)

#' @title Construct a Structured Quadratic Penalty
#'
#' @description
#' Builds the exact negative log-density of the Gaussian prior
#' \eqn{\beta \sim N(0, \Omega(\theta)^{-1})}, where the precision
#' \eqn{\Omega(\theta)} is a \pkg{parameters7} matrix parameter and the
#' hyperparameters are its free values. This is the correlated prior that
#' [quadratic_penalty()] cannot express: in a quadratic penalty one scale
#' multiplies a constant matrix, while here the hyperparameters reach every
#' entry of the matrix.
#'
#' @details
#' # The value
#'
#' Writing \eqn{\Omega} for the precision, which is the matrix that the
#' structure supplies,
#'
#' \deqn{\rho(\beta; \theta) = \tfrac{1}{2}\,\beta'\Omega(\theta)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,\Omega(\theta)
#'   + \tfrac{r}{2}\log 2\pi,}
#'
#' with \eqn{r} the structure's rank. For a full-rank structure this is the
#' negative log-density of the multivariate normal distribution with mean zero
#' and precision \eqn{\Omega(\theta)}; the examples compare it with a direct
#' evaluation of that density.
#'
#' # The structure is the precision
#'
#' The matrix that the structure supplies is the prior's precision:
#' \eqn{\Omega = M(\eta)}. This follows from the construction and is not an
#' option, because \eqn{\rho} is a negative log-density and the matrix that
#' enters its quadratic form is the precision.
#'
#' The choice of the side on which a structure is placed still matters,
#' because the two readings are different priors. For a positive definite
#' matrix they agree only where \eqn{M = M^{-1}}, that is at the identity. The
#' inverse of an AR(1) matrix is tridiagonal, so for a dimension of three or
#' more it is an AR(1) matrix only at zero correlation. To place a structure
#' `s` on the covariance side, build the family whose value is its inverse,
#'
#' ```
#' structured_penalty(parameters7::inverse_of(s))
#' ```
#'
#' which carries
#' \eqn{\partial_k \Sigma^{-1} = -\Sigma^{-1}B_k\Sigma^{-1}}, with
#' \eqn{B_k = \partial_k\Sigma} the derivative of the covariance, and its
#' higher orders inside the structure, where the log-determinant stays exact.
#' For [parameters7::ar1()] and [parameters7::autoregressive()] the inverse
#' has a constructor of its own, [parameters7::ar1_inv()] and
#' [parameters7::autoregressive_inv()], and a structure closed under inversion
#' requires no wrapper.
#'
#' # The derivatives
#'
#' Every derivative is computed from the derivative arrays that the structure
#' supplies. The hyperparameter gradient is
#' \eqn{\tfrac{1}{2}\beta'A_k\beta
#'   - \tfrac{1}{2}\partial_k\log\mathrm{pdet}\,\Omega},
#' with \eqn{A_k = \partial_k\Omega} the structure's `param_d1`; the Hessian
#' uses `param_d2`, and the mixed block is \eqn{A_k\beta}. The package does not
#' convert between covariance and precision: a structure meant for the
#' covariance side is a different structure, and it supplies its own arrays.
#'
#' # No map
#'
#' There is no `map` argument. A linear image of a structured precision is a
#' different precision, and composing it into the structure, where its
#' log-determinant stays exact, is left to \pkg{parameters7}.
#'
#' @param structure A \pkg{parameters7} `matrix_parameter`, read as the prior's
#'   precision. Pass `parameters7::inverse_of(s)` to put `s` on the covariance
#'   side instead.
#'
#'   It may be rank deficient, giving an improper prior: [is_proper()] then
#'   returns `FALSE` and the constant uses the rank and the log
#'   pseudo-determinant. A rank-deficient structure can be used only on the
#'   precision side, because `parameters7::inverse_of()` rejects it: a
#'   singular matrix has no inverse.
#'
#' @return A [StructuredPenalty()] object whose hyperparameters are the
#'   structure's `free_names`, each unconstrained and on an identity link.
#'
#' @examples
#' # An AR(1) prior on four coefficients: two hyperparameters reach every
#' # entry of the precision.
#' pen <- structured_penalty(parameters7::ar1(4))
#' pen@params
#' theta <- list(log_scale = 0.2, z_rho = 0.5)
#' b <- c(0.3, -0.1, 0.4, 0.2)
#' penalty_value(pen, b, theta)
#'
#' # The Hessian is the matrix of the structure, read as a precision.
#' Om <- parameters7::param_value(parameters7::ar1(4),
#'                                c(0.2, 0.5))
#' max(abs(penalty_hessian(pen, b, theta) - unclass(Om)))
#'
#' # The same structure read as a covariance is a different prior at the same
#' # free values: its Hessian is the inverse of that matrix.
#' cov <- structured_penalty(parameters7::inverse_of(parameters7::ar1(4)))
#' max(abs(penalty_hessian(cov, b, theta) - solve(unclass(Om))))
#' c(precision = penalty_value(pen, b, theta),
#'   covariance = penalty_value(cov, b, theta))
#'
#' # The covariance reading is the negative log-density of the multivariate
#' # normal with covariance S, written out directly.
#' S <- unclass(Om)
#' penalty_value(cov, b, theta) -
#'   (0.5 * sum(b * solve(S, b)) +
#'      0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) +
#'      2 * log(2 * pi))
#'
#' # At a zero log-Cholesky free vector the structure is the identity, so the
#' # penalty equals the plain ridge at lambda = 1, up to rounding.
#' s <- structured_penalty(parameters7::log_cholesky(3))
#' z <- as.list(stats::setNames(rep(0, length(s@params)), s@params))
#' bb <- c(0.4, -1.1, 0.7)
#' penalty_value(s, bb, z) - penalty_value(ridge_penalty(n_coef = 3), bb,
#'                                         list(lambda = 1))
#'
#' @seealso [quadratic_penalty()] for one scale on a fixed matrix,
#'   [additive_penalty()] for a sum of quadratics,
#'   [distrib_penalty()] for a coordinatewise prior,
#'   [parameters7::log_cholesky()] and [parameters7::ar1()] for structures to
#'   pass in, [penalty_readable()] for the interpretable quantities (a scale,
#'   a correlation) of a structure that declares them through
#'   [parameters7::param_readable()].
#' @export
structured_penalty <- function(structure) {
  if (!S7::S7_inherits(structure, parameters7::matrix_parameter)) {
    stop("'structure' must be a parameters7 matrix_parameter.", call. = FALSE)
  }
  nm <- structure@free_names
  StructuredPenalty(
    penalty_name = sprintf("structured [%s]", structure@param_name),
    map = NULL,
    n_coef = structure@dimension,
    params = nm,
    params_bounds = stats::setNames(rep(list(c(-Inf, Inf)), length(nm)), nm),
    link_params = stats::setNames(
      replicate(length(nm), linkfunctions7::identity_link(),
                simplify = FALSE), nm),
    params_smooth = stats::setNames(rep(TRUE, length(nm)), nm),
    structure = structure
  )
}

#' The Structure's Free Vector From the Aligned Hyperparameters
#'
#' @description
#' Unlists the aligned hyperparameter list into the numeric vector that the
#' generics of the structure take, in `pen@params` order, which is the order
#' of the structure's `free_names`.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param theta The aligned hyperparameter list, as [align_ptheta()] returns
#'   it.
#'
#' @return An unnamed numeric vector of length `length(pen@params)`.
#'
#' @seealso [struct_omega()], [align_ptheta()]
#'
#' @keywords internal
struct_eta <- function(pen, theta) {
  unlist(theta[pen@params], use.names = FALSE)
}

# The structure IS the precision: the helpers below unwrap its own answers and
# convert nothing. A covariance-side prior is a different structure,
# parameters7::inverse_of(s), which carries its own derivatives.
#
# The three matrix helpers return, respectively, Omega, the list of
# dOmega/dtheta_k in the order of pen@params, and the second derivatives keyed
# as param_d2 keys them ("name_i:name_j", the pair sorted by position and the
# key CONSTRUCTED from it, never parsed out of a name).

#' Precision, and Its Derivatives, From the Structure
#'
#' @description
#' This page documents the four helpers that read the prior's precision and its
#' derivatives off the structure. `struct_omega()` returns \eqn{\Omega},
#' `struct_d1()` its first derivatives, `struct_d2()` its second, and
#' `struct_logdet()` the log-determinant with its derivatives up to the
#' requested order.
#'
#' @details
#' The structure is the precision, so each helper removes the class from the
#' structure's own result, returns that result and does not convert between
#' covariance and precision. To read a structure `s` as a covariance, pass
#' `parameters7::inverse_of(s)`, which carries
#' \eqn{\partial_k\Omega = -\Omega B_k \Omega}, with \eqn{B_k} the derivative
#' of the covariance, and its higher orders inside the structure, where the
#' log-determinant stays exact.
#'
#' Second-order components are keyed as \pkg{parameters7} keys them: a key
#' joins the two free names with a colon, and the pair is sorted by position in
#' `free_names`. The key is constructed from the pair, never parsed out of a
#' name, since a free name may itself contain a colon.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param eta The structure's free vector, as [struct_eta()] returns it.
#' @param omega Accepted and unused; the calls in this package never pass it.
#'   Only `struct_d1()` and `struct_d2()` take it.
#' @param order The highest log-determinant derivative wanted: `0`, `1` or
#'   `2`. Only `struct_logdet()` takes it, and its default is `2L`.
#'
#' @return `struct_omega()` returns a symmetric base matrix of side
#'   `pen@structure@dimension`.
#'   `struct_d1()` returns a list of such matrices, one per free value, in
#'   `pen@params` order.
#'   `struct_d2()` returns a list of such matrices, one per unordered pair,
#'   keyed by the colon-joined sorted names. `struct_logdet()` returns a list
#'   with `value` (a single number) and, according to `order`, `d1` (a named
#'   numeric vector, one per free value) and `d2` (a named list, one number per
#'   pair, keyed as `struct_d2()` is).
#'
#' @seealso [structured_penalty()] for the prior that these serve,
#'   [parameters7::param_d1()] and [parameters7::param_dlogdet()] for the
#'   generics that they call.
#'
#' @keywords internal
struct_omega <- function(pen, eta) {
  unclass(parameters7::param_value(pen@structure, eta))
}

#' @rdname struct_omega
#' @keywords internal
struct_d1 <- function(pen, eta, omega = NULL) {
  lapply(parameters7::param_d1(pen@structure, eta), unclass)
}

#' @rdname struct_omega
#' @keywords internal
struct_d2 <- function(pen, eta, omega = NULL) {
  lapply(parameters7::param_d2(pen@structure, eta), unclass)
}

#' @rdname struct_omega
#' @keywords internal
struct_logdet <- function(pen, eta, order = 2L) {
  s <- pen@structure
  out <- list(value = parameters7::param_logdet(s, eta))
  if (order >= 1L) {
    out$d1 <- unlist(parameters7::param_dlogdet(s, eta))
  }
  if (order >= 2L) {
    out$d2 <- lapply(parameters7::param_d2logdet(s, eta), function(z) z)
  }
  out
}

#' @title Value of a Structured Penalty
#' @name penalty_value.StructuredPenalty
#'
#' @description
#' Returns the negative log-density of the Gaussian prior whose precision
#' \eqn{\Omega(\theta)} is the matrix of the structure. A structure built with
#' `parameters7::inverse_of()` supplies the inverse of the matrix that it
#' wraps, so that matrix is then the covariance.
#'
#' @details
#' \deqn{\rho(\beta; \theta) = \tfrac{1}{2}\,\beta'\Omega(\theta)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,\Omega(\theta)
#'   + \tfrac{r}{2}\log 2\pi,}
#'
#' with \eqn{r} the structure's rank. For a full-rank structure this is the
#' negative log-density of a multivariate normal distribution with mean zero
#' and precision \eqn{\Omega(\theta)}.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A named list of the structure's free values, already aligned by
#'   the generic. The bounds are the whole line, so this method rejects no
#'   value.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number.
#'
#' @examples
#' # A structure wrapped in inverse_of() places the AR(1) matrix S on the
#' # covariance side, and the value is the negative log-density of N(0, S).
#' s <- parameters7::inverse_of(parameters7::ar1(4))
#' pen <- structured_penalty(s)
#' b <- c(0.3, -0.1, 0.4, 0.2)
#' S <- unclass(parameters7::param_value(parameters7::ar1(4), c(0.2, 0.5)))
#' penalty_value(pen, b, list(log_scale = 0.2, z_rho = 0.5)) -
#'   (0.5 * sum(b * solve(S, b)) +
#'      0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) +
#'      2 * log(2 * pi))
#'
#' @seealso [structured_penalty()] for the construction,
#'   [penalty_gradient.StructuredPenalty()] for the coefficient derivatives,
#'   [struct_omega()] for the helpers that read the precision from the
#'   structure.
#' @keywords internal
S7::method(penalty_value, StructuredPenalty) <- function(pen, beta, theta, ...) {
  eta <- struct_eta(pen, theta)
  om <- struct_omega(pen, eta)
  sum(beta * as.numeric(om %*% beta)) / 2 -
    struct_logdet(pen, eta, 0L)$value / 2 +
    pen@structure@rank / 2 * log(2 * pi)
}

#' @title Coefficient Derivatives of a Structured Penalty
#' @name penalty_gradient.StructuredPenalty
#'
#' @description
#' `penalty_gradient()` returns \eqn{\Omega(\theta)\beta} and
#' `penalty_hessian()` returns \eqn{\Omega(\theta)} itself. The value is a
#' quadratic form in \eqn{\beta}, so the Hessian does not depend on the
#' coefficients and the gradient is the Hessian applied to them.
#'
#' @details
#' Both methods evaluate \eqn{\Omega(\theta)} through the structure's
#' `param_value()`. When the structure is `parameters7::inverse_of(s)`, that
#' evaluation includes the inversion of the matrix of `s`, and each generic
#' repeats it.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of the structure's free values.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` returns an unnamed numeric vector of length
#'   `pen@n_coef`; `penalty_hessian()` the symmetric precision matrix
#'   \eqn{\Omega(\theta)}, carrying the structure's own dimnames (`v1`,
#'   `v2`, ... for the \pkg{parameters7} primitives). The quadratic and
#'   additive branches return matrices with no dimnames, so a comparison
#'   across branches requires `unname()`.
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4))
#' th <- list(log_scale = 0.2, z_rho = 0.5)
#' b <- c(0.3, -0.1, 0.4, 0.2)
#'
#' # The gradient is the constant Hessian applied to the coefficients. The
#' # matrix carries the structure's dimnames and the vector does not.
#' all.equal(penalty_gradient(pen, b, th),
#'           unname(drop(penalty_hessian(pen, b, th) %*% b)))
#' dimnames(penalty_hessian(pen, b, th))[[1]]
#'
#' # And the Hessian is the structure's own matrix.
#' Om <- parameters7::param_value(parameters7::ar1(4),
#'                                c(0.2, 0.5))
#' max(abs(penalty_hessian(pen, b, th) - unclass(Om)))
#'
#' @seealso [penalty_value.StructuredPenalty()] for the quantity
#'   differentiated, [penalty_grad_theta.StructuredPenalty()] for the
#'   hyperparameter blocks.
#' @keywords internal
S7::method(penalty_gradient, StructuredPenalty) <- function(pen, beta, theta, ...) {
  as.numeric(struct_omega(pen, struct_eta(pen, theta)) %*% beta)
}

#' @rdname penalty_gradient.StructuredPenalty
#' @name penalty_hessian.StructuredPenalty
#' @keywords internal
S7::method(penalty_hessian, StructuredPenalty) <- function(pen, beta, theta, ...) {
  struct_omega(pen, struct_eta(pen, theta))
}

#' @title Hyperparameter Derivatives of a Structured Penalty
#' @name penalty_grad_theta.StructuredPenalty
#'
#' @description
#' The three blocks in the structure's free values are all assembled from the
#' structure's own derivative arrays. `penalty_grad_theta()` returns one number
#' per free value, `penalty_hess_theta()` one per unordered pair, and
#' `penalty_cross()` one coefficient vector per free value.
#'
#' @details
#' With \eqn{A_k = \partial_k\Omega} and \eqn{A_{kl} = \partial_{kl}\Omega} the
#' derivative arrays of the precision, as the structure supplies them,
#'
#' \deqn{\frac{\partial\rho}{\partial\theta_k}
#'     = \tfrac{1}{2}\beta'A_k\beta
#'       - \tfrac{1}{2}\partial_k \log\mathrm{pdet}\,\Omega, \qquad
#'   \frac{\partial^2\rho}{\partial\theta_k\partial\theta_l}
#'     = \tfrac{1}{2}\beta'A_{kl}\beta
#'       - \tfrac{1}{2}\partial_{kl} \log\mathrm{pdet}\,\Omega, \qquad
#'   \frac{\partial^2\rho}{\partial\beta\,\partial\theta_k} = A_k\beta.}
#'
#' The mixed block has no log-determinant term, because that term does not
#' depend on the coefficients.
#'
#' The second-order components are keyed by pairs as [penalty_hess_theta()]
#' keys them, diagonals first and then the upper off-diagonal pairs joined by
#' an underscore. Internally they are looked up by the structure's own keys,
#' each of which joins the two free names with a colon after sorting them by
#' position, and the key is constructed from the pair, not parsed out of a
#' name.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of the structure's free values.
#' @param scale Read by the generic, which applies the chain rule after
#'   dispatch. Every link here is the identity, so the two scales coincide and
#'   `"link"` returns the same numbers.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` returns a list of one number per free value,
#'   named by `pen@params`. `penalty_hess_theta()` returns a list of one number
#'   per unordered pair, keyed diagonals first. `penalty_cross()` returns a list
#'   of one numeric vector of length `pen@n_coef` per free value, named by
#'   `pen@params`.
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4))
#' th <- list(log_scale = 0.2, z_rho = 0.5)
#' b <- c(0.3, -0.1, 0.4, 0.2)
#'
#' penalty_grad_theta(pen, b, th)
#' names(penalty_hess_theta(pen, b, th))
#'
#' # Every link is the identity, so the unconstrained scale is the same one.
#' all.equal(penalty_grad_theta(pen, b, th),
#'           penalty_grad_theta(pen, b, th, scale = "link"))
#'
#' # The mixed block is A_k beta, so half its inner product with beta is the
#' # quadratic part of the gradient in z_rho, b' A_k b / 2.
#' cr <- penalty_cross(pen, b, th)
#' sum(b * cr$z_rho) / 2
#' sum(b * (parameters7::param_d1(parameters7::ar1(4), c(0.2, 0.5))$z_rho %*%
#'            b)) / 2
#'
#' @seealso [penalty_value.StructuredPenalty()] for the quantity
#'   differentiated, [struct_omega()] for the arrays these read,
#'   [penalty_logpdet.StructuredPenalty()] for the determinant's own
#'   derivatives.
#' @keywords internal
S7::method(penalty_grad_theta, StructuredPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    eta <- struct_eta(pen, theta)
    A <- struct_d1(pen, eta)
    dld <- struct_logdet(pen, eta, 1L)$d1
    stats::setNames(lapply(seq_along(pen@params), function(k) {
      sum(beta * as.numeric(A[[k]] %*% beta)) / 2 - dld[[k]] / 2
    }), pen@params)
  }

#' @rdname penalty_grad_theta.StructuredPenalty
#' @name penalty_hess_theta.StructuredPenalty
#' @keywords internal
S7::method(penalty_hess_theta, StructuredPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    eta <- struct_eta(pen, theta)
    A2 <- struct_d2(pen, eta)
    d2 <- struct_logdet(pen, eta, 2L)$d2
    nm <- pen@params
    prs <- ptheta_pairs(nm)
    # the structure keys its second-order components by the free names
    # joined with a colon, sorted by position; the key is CONSTRUCTED from
    # the pair, never parsed out of a name
    stats::setNames(lapply(prs, function(pr) {
      ij <- sort(match(pr, nm))
      key <- paste(nm[ij], collapse = ":")
      sum(beta * as.numeric(A2[[key]] %*% beta)) / 2 - d2[[key]] / 2
    }), names(prs))
  }

#' @rdname penalty_grad_theta.StructuredPenalty
#' @name penalty_cross.StructuredPenalty
#' @keywords internal
S7::method(penalty_cross, StructuredPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    A <- struct_d1(pen, struct_eta(pen, theta))
    stats::setNames(lapply(seq_along(pen@params), function(k) {
      as.numeric(A[[k]] %*% beta)
    }), pen@params)
  }

#' @title Smoothness and Kind of a Structured Penalty
#' @name penalty_kinks.StructuredPenalty
#'
#' @description
#' `penalty_kinks()` returns `numeric(0)`, because the value is a quadratic form
#' and smooth everywhere. `is_quadratic()` returns `TRUE`, so the four marginal
#' quantities are available. `is_proper()` returns `TRUE` when the structure
#' has full rank.
#'
#' @details
#' Properness is a comparison of the structure's rank with its dimension and
#' nothing else. A full-rank structure is a proper Gaussian prior. A
#' rank-deficient structure, which can be used only as a precision, gives an
#' improper prior, and its constant is the log pseudo-determinant over the
#' range.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param theta A named list of the structure's free values. It is
#'   accepted by `penalty_kinks()`, whose result does not depend on it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` returns a numeric vector of length zero.
#'   `is_quadratic()` returns the single logical `TRUE`.
#'   `is_proper()` returns a single logical, `TRUE` when the structure's rank
#'   equals its dimension.
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4))
#' penalty_kinks(pen, list(log_scale = 0.2, z_rho = 0.5))
#' is_quadratic(pen)
#' is_proper(pen)
#' c(penalty_rank(pen), pen@n_coef)
#'
#' @seealso [penalty_kinks()], [is_proper()] and [is_quadratic()] for the
#'   generics, [penalty_matrix.StructuredPenalty()] for the quantities that
#'   `is_quadratic()` gates.
#' @keywords internal
S7::method(penalty_kinks, StructuredPenalty) <- function(pen, theta, ...) {
  numeric(0)
}

#' @rdname penalty_kinks.StructuredPenalty
#' @name is_proper.StructuredPenalty
#' @keywords internal
S7::method(is_proper, StructuredPenalty) <- function(pen, ...) {
  pen@structure@rank == pen@structure@dimension
}

#' @rdname penalty_kinks.StructuredPenalty
#' @name is_quadratic.StructuredPenalty
#' @keywords internal
S7::method(is_quadratic, StructuredPenalty) <- function(pen, ...) TRUE

#' @title Marginal Quantities of a Structured Penalty
#' @name penalty_matrix.StructuredPenalty
#'
#' @description
#' This page documents the four quantities that a REML or marginal likelihood
#' criterion reads. `penalty_matrix()` returns the precision
#' \eqn{\Omega(\theta)}, `penalty_rank()` the structure's rank,
#' `penalty_null_basis()` its null basis, and `penalty_logpdet()` the log
#' pseudo-determinant of the precision with its first two derivatives in the
#' free values.
#'
#' @details
#' For a plain quadratic penalty \eqn{\log^{+}\lvert S\rvert} is linear in
#' \eqn{\log\lambda}. Here the log-determinant depends on every free value,
#' and its value and derivatives are those of the structure's
#' `param_logdet()`, `param_dlogdet()` and `param_d2logdet()`.
#'
#' The rank and the null basis are properties of the structure and do not
#' depend on \eqn{\theta}: a `matrix_parameter` records them at construction,
#' so a family whose rank is deficient is deficient in the same directions at
#' every free vector. `penalty_null_basis()` reshapes a vector-valued null
#' basis into a one-column matrix, so its return is always a matrix.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param theta A named list of the structure's free values. Only
#'   `penalty_matrix()` and `penalty_logpdet()` read it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_matrix()` returns the symmetric precision matrix, of side
#'   `pen@n_coef`, carrying the structure's own dimnames.
#'   `penalty_rank()` returns a single integer, the structure's rank.
#'   `penalty_null_basis()` returns a matrix with `pen@n_coef` rows and one
#'   column per null direction, with no columns for a full-rank structure.
#'   `penalty_logpdet()` returns a list of `value` (a number), `grad` (a list of
#'   one number per free value) and `hess` (a list of one number per unordered
#'   pair, keyed diagonals first).
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4))
#' th <- list(log_scale = 0.2, z_rho = 0.5)
#'
#' penalty_rank(pen)
#' dim(penalty_null_basis(pen))
#'
#' # The matrix is the precision itself, so it is the Hessian of the value.
#' all.equal(penalty_matrix(pen, th), penalty_hessian(pen, c(0, 0, 0, 0), th))
#'
#' # And its log pseudo-determinant is the ordinary log-determinant here,
#' # because the structure is full rank.
#' lp <- penalty_logpdet(pen, th)
#' lp$value - as.numeric(
#'   determinant(penalty_matrix(pen, th), logarithm = TRUE)$modulus)
#' names(lp$hess)
#'
#' @seealso [penalty_matrix()] and its three siblings for the generics,
#'   [struct_omega()] for where the derivatives come from,
#'   [penalty_matrix.QuadraticPenalty()] for the branch whose log
#'   pseudo-determinant is linear in \eqn{\log\lambda}.
#' @keywords internal
S7::method(penalty_matrix, StructuredPenalty) <- function(pen, theta, ...) {
  struct_omega(pen, struct_eta(pen, theta))
}

#' @rdname penalty_matrix.StructuredPenalty
#' @name penalty_rank.StructuredPenalty
#' @keywords internal
S7::method(penalty_rank, StructuredPenalty) <- function(pen, ...) {
  as.integer(pen@structure@rank)
}

#' @rdname penalty_matrix.StructuredPenalty
#' @name penalty_null_basis.StructuredPenalty
#' @keywords internal
S7::method(penalty_null_basis, StructuredPenalty) <- function(pen, ...) {
  nb <- pen@structure@null_basis
  if (is.null(dim(nb))) matrix(nb, nrow = pen@n_coef) else nb
}

#' @rdname penalty_matrix.StructuredPenalty
#' @name penalty_logpdet.StructuredPenalty
#' @keywords internal
S7::method(penalty_logpdet, StructuredPenalty) <- function(pen, theta, ...) {
  eta <- struct_eta(pen, theta)
  ld <- struct_logdet(pen, eta, 2L)
  nm <- pen@params
  prs <- ptheta_pairs(nm)
  list(
    value = ld$value,
    grad = stats::setNames(as.list(unname(ld$d1)), nm),
    hess = stats::setNames(lapply(prs, function(pr) {
      ij <- sort(match(pr, nm))
      unname(ld$d2[[paste(nm[ij], collapse = ":")]])
    }), names(prs))
  )
}

#' @title Interpretable Quantities of a Structured Penalty
#' @name penalty_readable.StructuredPenalty
#'
#' @description
#' Returns the interpretable quantities that the structure declares through
#' [parameters7::param_readable()], in the shape that
#' [penalty_readable.DistribPenalty()] returns, and `NULL` for a structure
#' that declares none.
#'
#' @details
#' The hyperparameters of a structured penalty are the free values of its
#' structure, such as the logarithm of a scale and the inverse hyperbolic
#' tangent of a correlation for [parameters7::ar1()].
#' [parameters7::param_readable()] maps them to the quantities that the
#' structure declares. Among the shipped structures, [parameters7::ar1()] and
#' [parameters7::compound_symmetry()] declare a scale and a correlation,
#' [parameters7::scaled_matrix()] its multiplier, and
#' [parameters7::autoregressive()] the scale, the partial autocorrelations
#' and the autoregressive coefficients. The other structures, among them
#' [parameters7::log_cholesky()], [parameters7::inverse_of()],
#' [parameters7::ar1_inv()] and the composition wrappers, declare nothing, and
#' the method then returns `NULL`.
#'
#' The quantities describe the matrix of the structure, and that matrix is
#' the prior's precision. For `structured_penalty(parameters7::ar1(4))` the
#' scale is the diagonal entry of the precision and the correlation is the
#' ratio of an entry next to the diagonal to the diagonal entry; they are not
#' the variance and the correlation of the coefficients. For this reason the
#' method prefixes every name that [parameters7::param_readable()] gives with
#' `precision_`, so that `scale` and `rho` are reported as `precision_scale`
#' and `precision_rho`. A structure placed on the covariance side through
#' [parameters7::inverse_of()] declares nothing, so the method returns `NULL`
#' for it.
#'
#' Every link of the branch is the identity, so the hyperparameters are the
#' free values themselves, and the Jacobian that
#' [parameters7::param_readable()] computes in closed form is the Jacobian
#' with respect to the hyperparameters. The method names its rows by the
#' prefixed names and its columns by `pen@params`, and repeats the
#' structure's label once per quantity to form `block`.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param theta A named list of the structure's free values.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `NULL` when the structure declares no quantities. Otherwise the
#'   function returns a list of four: `value`, a numeric vector of the
#'   quantities, in which each name is the structure's own name with the prefix
#'   `precision_`; `jacobian`, their derivatives, one row per quantity and one
#'   column per hyperparameter, with row names from `value` and column names
#'   from `pen@params`; `transform`, a character vector named like `value` that
#'   names the scale on which each interval is built; and `block`, the
#'   structure's label for every quantity, named like `value`.
#'
#' @examples
#' # The scale and the correlation of an AR(1) precision on four coefficients
#' # are named precision_scale and precision_rho.
#' pen <- structured_penalty(parameters7::ar1(4))
#' r <- penalty_readable(pen, list(log_scale = 0.2, z_rho = 0.5))
#' r$value
#' r$jacobian
#'
#' # The quantities describe the precision: its diagonal entry, and the ratio
#' # of an entry next to the diagonal to the diagonal entry.
#' Om <- penalty_matrix(pen, list(log_scale = 0.2, z_rho = 0.5))
#' c(Om[1, 1], Om[1, 2] / Om[1, 1])
#'
#' # A log-Cholesky structure declares nothing.
#' lc <- structured_penalty(parameters7::log_cholesky(2))
#' penalty_readable(lc, list(log_L1 = 0, log_L2 = 0, L2.1 = 0.3))
#'
#' @seealso [penalty_readable()] for the generic,
#'   [parameters7::param_readable()] for the quantities that a structure
#'   declares.
#' @keywords internal
S7::method(penalty_readable, StructuredPenalty) <- function(pen, theta, ...) {
  rd <- parameters7::param_readable(pen@structure, struct_eta(pen, theta))
  if (is.null(rd)) return(NULL)
  # the quantities describe the precision, not the coefficients' covariance
  nm <- paste0("precision_", names(rd$value))
  jac <- rd$jacobian
  dimnames(jac) <- list(nm, pen@params)
  list(value = stats::setNames(rd$value, nm), jacobian = jac,
       transform = stats::setNames(rd$transform, nm),
       block = stats::setNames(rep(rd$label, length(nm)), nm))
}
