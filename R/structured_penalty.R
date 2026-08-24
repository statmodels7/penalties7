#' @include quadratic_penalty.R
NULL

#' @title S7 Class for the Structured Quadratic Penalty
#'
#' @description
#' The class [structured_penalty()] builds: the Gaussian prior whose
#' covariance or precision is a \pkg{parameters7} matrix parameter. It adds one
#' property to [penalty()], the structure itself, and every quantity the branch
#' supplies is read off that structure's own contract.
#'
#' @details
#' The hyperparameters are the structure's free vector, so `params` is its
#' `free_names`, every entry of `params_bounds` is
#' \eqn{(-\infty, \infty)} and every link is the identity. The constraint that
#' keeps the matrix positive definite lives inside the structure, where a
#' scalar link cannot express it. `map` is always `NULL`.
#'
#' @inheritParams penalty
#' @param structure A \pkg{parameters7} `matrix_parameter` whose `role` is
#'   `"covariance"` or `"precision"`. It supplies the dimension, the rank, the
#'   null basis, the free names, the matrix and its derivative arrays.
#'
#' @return An S7 object of class `StructuredPenalty`, inheriting from
#'   [penalty()], with the seven inherited properties and `structure`.
#'
#' @seealso [structured_penalty()] for the constructor to use,
#'   [penalty_value.StructuredPenalty()] for what the branch computes,
#'   [quadratic_penalty()] for the branch whose matrix is fixed.
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
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
#' \eqn{\beta \sim N(0, \Sigma(\theta))}, where the matrix is a
#' \pkg{parameters7} matrix parameter and the hyperparameters are its free
#' values. This is the correlated prior [quadratic_penalty()] cannot express:
#' there one scale multiplies a constant matrix, here the hyperparameters reach
#' every entry.
#'
#' @details
#' # The value
#'
#' Writing \eqn{\Omega} for the precision, whether the structure supplies it
#' directly or through the covariance it describes,
#'
#' \deqn{\rho(\beta; \theta) = \tfrac{1}{2}\,\beta'\Omega(\theta)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,\Omega(\theta)
#'   + \tfrac{r}{2}\log 2\pi,}
#'
#' with \eqn{r} the structure's rank. For a full-rank structure this is exactly
#' `-mvtnorm`-style multivariate normal log-density, which the examples check
#' against a hand-written one.
#'
#' # The two roles
#'
#' The structure's `role` says which matrix of the prior it is, and the two
#' readings are different priors from the same free vector. A structure of
#' role `"precision"` at \eqn{\eta} gives \eqn{\Omega = M(\eta)}; one of role
#' `"covariance"` gives \eqn{\Omega = M(\eta)^{-1}}. There is in general no
#' free vector that makes the two agree: the inverse of an AR(1) covariance is
#' tridiagonal and is not an AR(1) covariance at any parameters.
#'
#' A structure that declares `"either"` is rejected. The two readings differ in
#' the sign of the log-determinant term, and nothing in the matrix says which
#' was meant, so a default would give a fit that converges to a different prior
#' without saying so.
#'
#' # The derivatives
#'
#' Every derivative comes from the structure's own contract. Where the
#' structure is the precision, the hyperparameter gradient is
#' \eqn{\tfrac{1}{2}\beta'A_k\beta - \tfrac{1}{2}\partial_k\log\mathrm{pdet}}
#' with \eqn{A_k} the structure's `param_d1`, the Hessian adds `param_d2`, and
#' the mixed block is \eqn{A_k\beta}. Where it is the covariance the same
#' expressions are read at the precision it implies, whose derivatives follow
#' from the chain rule for an inverse,
#'
#' \deqn{\partial_k\Omega = -\Omega A_k \Omega, \qquad
#'   \partial_{kl}\Omega = \Omega\left(A_k\Omega A_l
#'     + A_l\Omega A_k\right)\Omega - \Omega A_{kl}\Omega,}
#'
#' with \eqn{\log\lvert\Omega\rvert = -\log\lvert\Sigma\rvert} and its
#' derivatives negated termwise. The transport is done once, in
#' [struct_omega()] and its siblings, so the methods are the same arithmetic in
#' both cases.
#'
#' # No map
#'
#' There is no `map` argument. A linear image of a structured precision is a
#' different precision, and composing it into the structure, where its
#' log-determinant stays exact, is the structure's own business.
#'
#' @param structure A \pkg{parameters7} `matrix_parameter` whose `role` says
#'   which matrix of the prior it is.
#'
#'   A structure declared `"precision"` may be rank deficient, giving an
#'   improper prior: [is_proper()] then answers `FALSE` and the constant uses
#'   the rank and the log pseudo-determinant.
#'
#'   A structure declared `"covariance"` may not be rank deficient. A
#'   covariance of deficient rank has no inverse, and a direction of zero
#'   variance is a constraint on the coefficients and not a prior over them;
#'   the constructor rejects it, naming the rank and the dimension.
#'
#'   A structure declaring `"either"` is rejected, with the two roles named.
#'
#' @return A [StructuredPenalty()] object whose hyperparameters are the
#'   structure's `free_names`, each unconstrained and on an identity link.
#'
#' @examples
#' # An AR(1) prior on four coefficients: two hyperparameters reach every
#' # entry of the precision.
#' pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
#' pen@params
#' theta <- list(log_scale = 0.2, z_rho = 0.5)
#' b <- c(0.3, -0.1, 0.4, 0.2)
#' penalty_value(pen, b, theta)
#'
#' # The Hessian is the structure's own matrix, read as a precision.
#' Om <- parameters7::param_value(parameters7::ar1(4, role = "precision"),
#'                                c(0.2, 0.5))
#' max(abs(penalty_hessian(pen, b, theta) - unclass(Om)))
#'
#' # The same structure read as a covariance is a DIFFERENT prior at the same
#' # free values: there the Hessian is the inverse of that matrix.
#' cov <- structured_penalty(parameters7::ar1(4, role = "covariance"))
#' max(abs(penalty_hessian(cov, b, theta) - solve(unclass(Om))))
#' c(precision = penalty_value(pen, b, theta),
#'   covariance = penalty_value(cov, b, theta))
#'
#' # And the covariance reading is exactly a multivariate normal log-density.
#' S <- unclass(Om)
#' penalty_value(cov, b, theta) -
#'   (0.5 * sum(b * solve(S, b)) +
#'      0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) +
#'      2 * log(2 * pi))
#'
#' # At a zero log-Cholesky free vector the structure is the identity, so the
#' # penalty is the plain ridge at lambda = 1, to the last bit.
#' s <- structured_penalty(parameters7::log_cholesky(3, role = "precision"))
#' z <- as.list(stats::setNames(rep(0, length(s@params)), s@params))
#' bb <- c(0.4, -1.1, 0.7)
#' penalty_value(s, bb, z) - penalty_value(ridge_penalty(n_coef = 3), bb,
#'                                         list(lambda = 1))
#'
#' @seealso [quadratic_penalty()] for one scale on a fixed matrix,
#'   [additive_penalty()] for a sum of quadratics,
#'   [distrib_penalty()] for a coordinatewise prior,
#'   [parameters7::log_cholesky()] and [parameters7::ar1()] for structures to
#'   pass in, [penalty_readable()] for reporting the hyperparameters as
#'   standard deviations and correlations.
#' @export
structured_penalty <- function(structure) {
  if (!S7::S7_inherits(structure, parameters7::matrix_parameter)) {
    stop("'structure' must be a parameters7 matrix_parameter.", call. = FALSE)
  }
  # The role is READ, not defaulted. A structure that serves as either is a
  # statement about the structure and not about this prior, and the two
  # readings differ in the sign of the log-determinant term: guessing would
  # give a fit that converges to a different matrix without saying so.
  role <- structure@role
  if (!identical(role, "covariance") && !identical(role, "precision")) {
    stop(sprintf(paste0(
      "'%s' declares role '%s', so which matrix of the prior it is has not\n",
      "  been said. Rebuild it with role = \"covariance\" or\n",
      "  role = \"precision\": the two differ in the sign of the\n",
      "  log-determinant term and cannot be told apart from the matrix."),
      structure@param_name, role), call. = FALSE)
  }
  if (identical(role, "covariance") && structure@rank < structure@dimension) {
    stop(sprintf(paste0(
      "'%s' is a covariance of rank %d out of %d, so it has no inverse and\n",
      "  the prior does not exist: a direction of zero variance is a\n",
      "  constraint on the coefficients, not a prior over them. A\n",
      "  rank-deficient structure is admitted as a PRECISION, where it is\n",
      "  the improper prior the log pseudo-determinant is written for."),
      structure@param_name, structure@rank, structure@dimension),
      call. = FALSE)
  }
  nm <- structure@free_names
  StructuredPenalty(
    penalty_name = sprintf("structured [%s, %s]", structure@param_name, role),
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

#' Whether the Structure Describes the Covariance
#'
#' @description
#' Reads the structure's declared role and answers `TRUE` for
#' `"covariance"`. Every method of the branch is written in the precision, so
#' this is the one place that decides whether a transport is needed.
#'
#' @param pen A [StructuredPenalty()] object.
#'
#' @return A single logical. `FALSE` for a structure of role `"precision"`,
#'   which is the only other value [structured_penalty()] admits.
#'
#' @seealso [struct_omega()] for the transport this gates
#'
#' @keywords internal
struct_is_cov <- function(pen) identical(pen@structure@role, "covariance")

#' The Structure's Free Vector From the Aligned Hyperparameters
#'
#' @description
#' Unlists the aligned hyperparameter list into the numeric vector the
#' structure's own generics take, in `pen@params` order, which is the
#' structure's `free_names` order.
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

# The penalty is written in the PRECISION whatever the structure describes, so
# the transport happens here and every method below is the same arithmetic in
# both cases. Writing it twice would be two implementations of one prior, and
# the second one is the one nobody reads.
#
# The three helpers return, respectively, Omega, the list of dOmega/dtheta_k in
# the order of pen@params, and the second derivatives keyed as param_d2 keys
# them ("name_i:name_j", the pair sorted by position and the key CONSTRUCTED
# from it, never parsed out of a name).

#' Precision, and Its Derivatives, From the Structure
#'
#' @description
#' The four helpers that read the prior's precision and its derivatives off the
#' structure, transporting from the covariance where that is what the structure
#' describes. `struct_omega()` returns \eqn{\Omega}, `struct_d1()` its first
#' derivatives, `struct_d2()` its second, and `struct_logdet()` the
#' log-determinant with as many orders as asked for.
#'
#' @details
#' The branch is written in the precision throughout, so the transport happens
#' here and every method is the same arithmetic in both roles. Writing it twice
#' would be two implementations of one prior.
#'
#' For a covariance structure, with \eqn{A_k} and \eqn{A_{kl}} the structure's
#' own derivative arrays,
#'
#' \deqn{\Omega = \Sigma^{-}, \qquad \partial_k\Omega = -\Omega A_k \Omega,
#'   \qquad \partial_{kl}\Omega = \Omega\left(A_k\Omega A_l
#'     + A_l\Omega A_k\right)\Omega - \Omega A_{kl}\Omega,}
#'
#' and the log-determinant is negated at every order. For a precision structure
#' each helper unwraps the structure's answer and returns it.
#'
#' Second-order components are keyed as \pkg{parameters7} keys them: the two
#' free names joined by a colon, the pair sorted by position in `free_names`.
#' The key is constructed from the pair, never parsed out of a name, since a
#' free name may itself contain a colon.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param eta The structure's free vector, as [struct_eta()] returns it.
#' @param omega The precision, when the caller already has it, so that a
#'   covariance structure is not inverted twice. `NULL` computes it.
#'   `struct_d1()` and `struct_d2()` only.
#' @param order The highest log-determinant derivative wanted: `0`, `1` or
#'   `2`. `struct_logdet()` only, `2L` by default.
#'
#' @return `struct_omega()` a symmetric base matrix of side
#'   `pen@structure@dimension`.
#'   `struct_d1()` a list of such matrices, one per free value, in
#'   `pen@params` order.
#'   `struct_d2()` a list of such matrices, one per unordered pair, keyed by
#'   the colon-joined sorted names.
#'   `struct_logdet()` a list with `value` (a single number) and, according to
#'   `order`, `d1` (a named numeric vector, one per free value) and `d2` (a
#'   named list, one number per pair, keyed as `struct_d2()` is).
#'
#' @seealso [structured_penalty()] for the prior these serve,
#'   [parameters7::param_d1()] and [parameters7::param_dlogdet()] for the
#'   contract they read.
#'
#' @keywords internal
struct_omega <- function(pen, eta) {
  s <- pen@structure
  if (struct_is_cov(pen)) {
    unclass(parameters7::param_solve(s, eta))
  } else {
    unclass(parameters7::param_value(s, eta))
  }
}

#' @rdname struct_omega
#' @keywords internal
struct_d1 <- function(pen, eta, omega = NULL) {
  A <- parameters7::param_d1(pen@structure, eta)
  if (!struct_is_cov(pen)) return(lapply(A, unclass))
  if (is.null(omega)) omega <- struct_omega(pen, eta)
  lapply(A, function(Ak) -(omega %*% unclass(Ak) %*% omega))
}

#' @rdname struct_omega
#' @keywords internal
struct_d2 <- function(pen, eta, omega = NULL) {
  s <- pen@structure
  A2 <- parameters7::param_d2(s, eta)
  if (!struct_is_cov(pen)) return(lapply(A2, unclass))
  if (is.null(omega)) omega <- struct_omega(pen, eta)
  A <- lapply(parameters7::param_d1(s, eta), unclass)
  nm <- pen@params
  prs <- ptheta_pairs(nm)
  out <- lapply(prs, function(pr) {
    ij <- sort(match(pr, nm))
    key <- paste(nm[ij], collapse = ":")
    Ak <- A[[ij[1L]]]
    Al <- A[[ij[2L]]]
    omega %*% (Ak %*% omega %*% Al + Al %*% omega %*% Ak) %*% omega -
      omega %*% unclass(A2[[key]]) %*% omega
  })
  names(out) <- vapply(prs, function(pr) {
    ij <- sort(match(pr, nm))
    paste(nm[ij], collapse = ":")
  }, "")
  out
}

#' @rdname struct_omega
#' @param order The highest derivative wanted, 0, 1 or 2.
#' @keywords internal
struct_logdet <- function(pen, eta, order = 2L) {
  s <- pen@structure
  sgn <- if (struct_is_cov(pen)) -1 else 1
  out <- list(value = sgn * parameters7::param_logdet(s, eta))
  if (order >= 1L) {
    out$d1 <- sgn * unlist(parameters7::param_dlogdet(s, eta))
  }
  if (order >= 2L) {
    out$d2 <- lapply(parameters7::param_d2logdet(s, eta), function(z) sgn * z)
  }
  out
}

#' @title Value of a Structured Penalty
#' @name penalty_value.StructuredPenalty
#'
#' @description
#' Returns the negative log-density of the Gaussian prior the structure
#' describes, evaluated in the precision whichever of the two matrices the
#' structure supplies.
#'
#' @details
#' \deqn{\rho(\beta; \theta) = \tfrac{1}{2}\,\beta'\Omega(\theta)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,\Omega(\theta)
#'   + \tfrac{r}{2}\log 2\pi,}
#'
#' with \eqn{r} the structure's rank, and the log-determinant negated where the
#' structure is a covariance. For a full-rank structure this is exactly a
#' multivariate normal log-density with mean zero.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A named list of the structure's free values, already aligned by
#'   the generic. The bounds are the whole line, so nothing is rejected here.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number.
#'
#' @examples
#' # A covariance structure gives exactly a multivariate normal log-density.
#' s <- parameters7::ar1(4, role = "covariance")
#' pen <- structured_penalty(s)
#' b <- c(0.3, -0.1, 0.4, 0.2)
#' S <- unclass(parameters7::param_value(s, c(0.2, 0.5)))
#' penalty_value(pen, b, list(log_scale = 0.2, z_rho = 0.5)) -
#'   (0.5 * sum(b * solve(S, b)) +
#'      0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus) +
#'      2 * log(2 * pi))
#'
#' @seealso [structured_penalty()] for the construction,
#'   [penalty_gradient.StructuredPenalty()] for the coefficient derivatives,
#'   [struct_omega()] for the transport.
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
#' quadratic form in \eqn{\beta}, so the Hessian does not move with the
#' coefficients and the gradient is the Hessian applied to them.
#'
#' @details
#' Where the structure is a covariance, both read the inverse it implies, so
#' each call costs one solve. A caller taking several derivatives at one
#' \eqn{\theta} pays for that inversion once per generic.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of the structure's free values.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` an unnamed numeric vector of length
#'   `pen@n_coef`; `penalty_hessian()` the symmetric precision matrix of that
#'   side, **carrying the structure's own dimnames** (`v1`, `v2`, ... for the
#'   \pkg{parameters7} primitives). The other branches return matrices with no
#'   dimnames, so a comparison across branches wants `unname()`.
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
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
#' Om <- parameters7::param_value(parameters7::ar1(4, role = "precision"),
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
#' The three blocks in the structure's free values, all assembled from the
#' structure's own derivative arrays. `penalty_grad_theta()` returns one number
#' per free value, `penalty_hess_theta()` one per unordered pair, and
#' `penalty_cross()` one coefficient vector per free value.
#'
#' @details
#' With \eqn{A_k = \partial_k\Omega} and \eqn{A_{kl} = \partial_{kl}\Omega} the
#' precision's derivative arrays, transported from the covariance where the
#' structure describes one,
#'
#' \deqn{\frac{\partial\rho}{\partial\theta_k}
#'     = \tfrac{1}{2}\beta'A_k\beta
#'       - \tfrac{1}{2}\partial_k \log\mathrm{pdet}\,\Omega, \qquad
#'   \frac{\partial^2\rho}{\partial\theta_k\partial\theta_l}
#'     = \tfrac{1}{2}\beta'A_{kl}\beta
#'       - \tfrac{1}{2}\partial_{kl} \log\mathrm{pdet}\,\Omega, \qquad
#'   \frac{\partial^2\rho}{\partial\beta\,\partial\theta_k} = A_k\beta.}
#'
#' The mixed block carries no log-determinant term, that term not depending on
#' the coefficients.
#'
#' The second-order components are keyed by pairs as [penalty_hess_theta()]
#' keys them, diagonals first and then the upper off-diagonal pairs joined by
#' an underscore. Internally they are looked up in the structure's own keying,
#' the two free names joined by a colon and sorted by position, and the key is
#' constructed from the pair rather than parsed out of a name.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of the structure's free values.
#' @param scale Read by the generic, which applies the chain rule after
#'   dispatch. Every link here is the identity, so the two scales coincide and
#'   `"link"` returns the same numbers.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` a list of one number per free value, named by
#'   `pen@params`.
#'   `penalty_hess_theta()` a list of one number per unordered pair, keyed
#'   diagonals first.
#'   `penalty_cross()` a list of one numeric vector of length `pen@n_coef` per
#'   free value, named by `pen@params`.
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
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
#' # The mixed block is A_k beta, so summing it against beta gives twice the
#' # quadratic part of the gradient in that direction.
#' cr <- penalty_cross(pen, b, th)
#' sum(b * cr$z_rho) / 2
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
#' `penalty_kinks()` returns `numeric(0)`, the value being a quadratic form and
#' smooth everywhere. `is_quadratic()` returns `TRUE`, so the four marginal
#' quantities are available. `is_proper()` returns `TRUE` when the structure
#' has full rank.
#'
#' @details
#' Properness is the structure's rank against its dimension and nothing else. A
#' full-rank structure is a proper Gaussian prior; a rank-deficient one, which
#' is admitted only as a precision, is improper and its constant is the log
#' pseudo-determinant over the range.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param theta A named list of the structure's free values. Read by
#'   `penalty_kinks()`, whose answer does not depend on it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` a numeric vector of length zero.
#'   `is_quadratic()` the single logical `TRUE`.
#'   `is_proper()` a single logical, `TRUE` when the structure's rank equals
#'   its dimension.
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
#' penalty_kinks(pen, list(log_scale = 0.2, z_rho = 0.5))
#' is_quadratic(pen)
#' is_proper(pen)
#' c(penalty_rank(pen), pen@n_coef)
#'
#' @seealso [penalty_kinks()], [is_proper()] and [is_quadratic()] for the
#'   generics, [penalty_matrix.StructuredPenalty()] for the quantities
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
#' The four pieces a REML or marginal likelihood criterion reads.
#' `penalty_matrix()` returns the precision \eqn{\Omega(\theta)},
#' `penalty_rank()` the structure's rank, `penalty_null_basis()` its null
#' basis, and `penalty_logpdet()` the log pseudo-determinant of the precision
#' with its first two derivatives in the free values.
#'
#' @details
#' Where a plain quadratic penalty has \eqn{\log^{+}\lvert S\rvert} linear in
#' \eqn{\log\lambda}, here the log-determinant moves with every free value and
#' its derivatives come from the structure's own log-determinant contract,
#' negated at every order when the structure describes a covariance.
#'
#' The rank and the null basis are properties of the structure and do not
#' depend on \eqn{\theta}: a `matrix_parameter` records them at construction,
#' so a family whose rank is deficient is deficient in the same directions at
#' every free vector. `penalty_null_basis()` reshapes a vector-valued null
#' basis into a one-column matrix, so its return is always a matrix.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param theta A named list of the structure's free values.
#'   `penalty_matrix()` and `penalty_logpdet()` only.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_matrix()` the symmetric precision matrix, of side
#'   `pen@n_coef`, carrying the structure's own dimnames.
#'   `penalty_rank()` a single integer, the structure's rank.
#'   `penalty_null_basis()` a matrix with `pen@n_coef` rows and one column per
#'   null direction, with no columns for a full-rank structure.
#'   `penalty_logpdet()` a list of `value` (a number), `grad` (a list of one
#'   number per free value) and `hess` (a list of one number per unordered
#'   pair, keyed diagonals first).
#'
#' @examples
#' pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
#' th <- list(log_scale = 0.2, z_rho = 0.5)
#'
#' penalty_rank(pen)
#' dim(penalty_null_basis(pen))
#'
#' # The matrix is the precision itself, so it is the Hessian of the value.
#' all.equal(penalty_matrix(pen, th), penalty_hessian(pen, c(0, 0, 0, 0), th))
#'
#' # And its log pseudo-determinant is the ordinary log-determinant here,
#' # the structure being full rank.
#' lp <- penalty_logpdet(pen, th)
#' lp$value - as.numeric(
#'   determinant(penalty_matrix(pen, th), logarithm = TRUE)$modulus)
#' names(lp$hess)
#'
#' @seealso [penalty_matrix()] and its three siblings for the generics,
#'   [struct_omega()] for where the derivatives come from,
#'   [penalty_matrix.QuadraticPenalty()] for the branch whose determinant is
#'   linear in one parameter.
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
