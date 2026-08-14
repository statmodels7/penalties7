#' @include quadratic_penalty.R
NULL

#' @title S7 Class for the Structured Quadratic Penalty
#'
#' @description
#' The class \code{\link{structured_penalty}} instantiates: the Gaussian
#' prior whose covariance or precision is a \pkg{parameters7} matrix
#' parameter, so that the hyperparameters enter the matrix itself.
#'
#' @inheritParams penalty
#' @param structure The \pkg{parameters7} \code{matrix_parameter}.
#'
#' @return An object of class \code{StructuredPenalty}.
#'
#' @seealso \code{\link{structured_penalty}}
#' @examples
#' S7::S7_inherits(
#'   structured_penalty(parameters7::log_cholesky(2, role = "precision")),
#'   StructuredPenalty)
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
#' The exact negative log-density of the Gaussian prior
#' \eqn{\beta \sim N(0, \Sigma(\theta))}, with the matrix a
#' \pkg{parameters7} matrix parameter. Writing \eqn{\Omega = \Sigma^{-}} for
#' the precision,
#' \deqn{\rho = \tfrac{1}{2}\,\beta'\Omega(\theta)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,\Omega(\theta)
#'   + \tfrac{r}{2}\log 2\pi,}
#' and the structure supplies whichever of the two matrices its \code{role}
#' declares. This is the correlated prior \code{\link{quadratic_penalty}}
#' cannot express: there the matrix is a constant and one scale multiplies
#' it, here the hyperparameters are the structure's free values and reach
#' every entry.
#'
#' @details
#' The hyperparameters ARE the structure's free vector, which is
#' unconstrained by construction, so every link is the identity -- the same
#' flattening convention the multivariate families of \pkg{distributions7}
#' follow, and for the same reason: the constraint lives in the structure,
#' where a scalar link cannot express it.
#'
#' Every derivative comes from the structure's own contract. When the
#' structure IS the precision the theta gradient is
#' \eqn{\tfrac{1}{2}\beta'A_k\beta - \tfrac{1}{2}\partial_k\log\mathrm{pdet}}
#' with \eqn{A_k} the structure's \code{param_d1}, the theta Hessian adds
#' \code{param_d2}, and the mixed block is \eqn{A_k\beta}. When it is the
#' covariance the same expressions are read at the precision it implies,
#' whose derivatives follow from the chain rule for an inverse,
#' \deqn{\partial_k\Omega = -\Omega A_k \Omega, \qquad
#'   \partial_{kl}\Omega = \Omega\left(A_k\Omega A_l
#'     + A_l\Omega A_k\right)\Omega - \Omega A_{kl}\Omega,}
#' with \eqn{\log\lvert\Omega\rvert = -\log\lvert\Sigma\rvert} and its
#' derivatives negated termwise. The transport is done once and the
#' quantities below are then the same arithmetic in both cases.
#'
#' There is no \code{map} argument, deliberately: a linear image of a
#' structured precision is a different precision, and composing it into the
#' structure -- where its log-determinant stays exact -- is the structure's
#' business, not this constructor's.
#'
#' @param structure A \pkg{parameters7} \code{matrix_parameter} whose
#'   \code{role} says which matrix of the prior it is. A structure declared
#'   \code{"precision"} may be rank deficient (an improper prior);
#'   \code{is_proper} then answers \code{FALSE} and the constant uses the
#'   rank and the log pseudo-determinant. A structure declared
#'   \code{"covariance"} may not: a covariance of deficient rank has no
#'   inverse, and an effect with a direction of zero variance is a
#'   constraint rather than a prior. A structure that declares
#'   \code{"either"} is rejected, because the sign of the log-determinant
#'   term depends on which of the two it is.
#'
#' @return An object of class \code{StructuredPenalty}.
#'
#' @examples
#' # an AR(1) prior on four coefficients: three hyperparameters reach every
#' # entry of the precision
#' pen <- structured_penalty(parameters7::ar1(4, role = "precision"))
#' theta <- list(log_scale = 0.2, z_rho = 0.5)
#' penalty_value(pen, c(0.3, -0.1, 0.4, 0.2), theta)
#' penalty_rank(pen)
#'
#' # the same prior written by its covariance: the free values that give
#' # Sigma here give Sigma^-1 there, and the penalty is the same number
#' cov <- structured_penalty(parameters7::ar1(4, role = "covariance"))
#' penalty_value(cov, c(0.3, -0.1, 0.4, 0.2), theta)
#'
#' @seealso \code{\link{quadratic_penalty}}, \code{\link{additive_penalty}}, \code{\link{distrib_penalty}}
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
#' @description Reads the structure's declared role.
#' @param pen A \code{StructuredPenalty} object.
#' @return A single logical.
#' @keywords internal
struct_is_cov <- function(pen) identical(pen@structure@role, "covariance")

#' The Structure's Free Vector From the Aligned Hyperparameters
#' @description Unlists the aligned theta in the structure's own order.
#' @param pen A \code{StructuredPenalty} object.
#' @param theta The aligned hyperparameter list.
#' @return A numeric vector.
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
#' @description
#' The prior's precision and its first and second derivatives in the
#' structure's free values, transported from the covariance where that is what
#' the structure describes.
#' @param pen A \code{StructuredPenalty} object.
#' @param eta The structure's free vector.
#' @param omega The precision, when the caller already has it.
#' @return A matrix (\code{struct_omega}) or a list of matrices.
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

#' @title Structured Penalty Methods
#' @name penalty_value.StructuredPenalty
#' @description
#' Every quantity from the structure's own contract; see
#' \code{\link{structured_penalty}}.
#' @param pen A \code{StructuredPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of the structure's free values.
#' @param scale Handled by the generic.
#' @param ... Unused.
#' @return See the generic pages.
#' @keywords internal
S7::method(penalty_value, StructuredPenalty) <- function(pen, beta, theta, ...) {
  eta <- struct_eta(pen, theta)
  om <- struct_omega(pen, eta)
  sum(beta * as.numeric(om %*% beta)) / 2 -
    struct_logdet(pen, eta, 0L)$value / 2 +
    pen@structure@rank / 2 * log(2 * pi)
}

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_gradient.StructuredPenalty
#' @keywords internal
S7::method(penalty_gradient, StructuredPenalty) <- function(pen, beta, theta, ...) {
  as.numeric(struct_omega(pen, struct_eta(pen, theta)) %*% beta)
}

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_hessian.StructuredPenalty
#' @keywords internal
S7::method(penalty_hessian, StructuredPenalty) <- function(pen, beta, theta, ...) {
  struct_omega(pen, struct_eta(pen, theta))
}

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_grad_theta.StructuredPenalty
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

#' @rdname penalty_value.StructuredPenalty
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

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_cross.StructuredPenalty
#' @keywords internal
S7::method(penalty_cross, StructuredPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    A <- struct_d1(pen, struct_eta(pen, theta))
    stats::setNames(lapply(seq_along(pen@params), function(k) {
      as.numeric(A[[k]] %*% beta)
    }), pen@params)
  }

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_kinks.StructuredPenalty
#' @keywords internal
S7::method(penalty_kinks, StructuredPenalty) <- function(pen, theta, ...) {
  numeric(0)
}

#' @rdname penalty_value.StructuredPenalty
#' @name is_proper.StructuredPenalty
#' @keywords internal
S7::method(is_proper, StructuredPenalty) <- function(pen, ...) {
  pen@structure@rank == pen@structure@dimension
}

#' @rdname penalty_value.StructuredPenalty
#' @name is_quadratic.StructuredPenalty
#' @keywords internal
S7::method(is_quadratic, StructuredPenalty) <- function(pen, ...) TRUE

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_matrix.StructuredPenalty
#' @keywords internal
S7::method(penalty_matrix, StructuredPenalty) <- function(pen, theta, ...) {
  struct_omega(pen, struct_eta(pen, theta))
}

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_rank.StructuredPenalty
#' @keywords internal
S7::method(penalty_rank, StructuredPenalty) <- function(pen, ...) {
  as.integer(pen@structure@rank)
}

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_null_basis.StructuredPenalty
#' @keywords internal
S7::method(penalty_null_basis, StructuredPenalty) <- function(pen, ...) {
  nb <- pen@structure@null_basis
  if (is.null(dim(nb))) matrix(nb, nrow = pen@n_coef) else nb
}

#' @rdname penalty_value.StructuredPenalty
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
