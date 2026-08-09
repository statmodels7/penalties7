#' @include quadratic_penalty.R
NULL

#' @title S7 Class for the Structured Quadratic Penalty
#'
#' @description
#' The class \code{\link{structured_penalty}} instantiates: the Gaussian
#' prior whose precision is a \pkg{parameters7} matrix parameter, so that
#' the hyperparameters enter the matrix itself.
#'
#' @inheritParams penalty
#' @param structure The \pkg{parameters7} \code{matrix_parameter}.
#'
#' @return An object of class \code{StructuredPenalty}.
#'
#' @seealso \code{\link{structured_penalty}}
#' @examples
#' S7::S7_inherits(structured_penalty(parameters7::log_cholesky(2)),
#'                 StructuredPenalty)
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
#' \eqn{\beta \sim N(0, \Omega(\theta)^{-})}, with the precision
#' \eqn{\Omega} a \pkg{parameters7} matrix parameter:
#' \deqn{\rho = \tfrac{1}{2}\,\beta'\Omega(\theta)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,\Omega(\theta)
#'   + \tfrac{r}{2}\log 2\pi.}
#' This is the correlated prior \code{\link{quadratic_penalty}} cannot
#' express: there the matrix is a constant and one scale multiplies it,
#' here the hyperparameters are the structure's free values and reach every
#' entry.
#'
#' @details
#' The hyperparameters ARE the structure's free vector, which is
#' unconstrained by construction, so every link is the identity -- the same
#' flattening convention the multivariate families of \pkg{distributions7}
#' follow, and for the same reason: the constraint lives in the structure,
#' where a scalar link cannot express it. Every derivative comes from the
#' structure's own contract: the theta gradient is
#' \eqn{\tfrac{1}{2}\beta'A_k\beta - \tfrac{1}{2}\partial_k\log\mathrm{pdet}}
#' with \eqn{A_k} the structure's \code{param_d1}, the theta Hessian adds
#' \code{param_d2}, and the mixed block is \eqn{A_k\beta}.
#'
#' There is no \code{map} argument, deliberately: a linear image of a
#' structured precision is a different precision, and composing it into the
#' structure -- where its log-determinant stays exact -- is the structure's
#' business, not this constructor's.
#'
#' @param structure A \pkg{parameters7} \code{matrix_parameter}, read as the
#'   PRECISION of the prior. A rank-deficient structure (an improper prior)
#'   is admitted; \code{is_proper} then answers \code{FALSE} and the
#'   constant uses the rank and the log pseudo-determinant.
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
#' @seealso \code{\link{quadratic_penalty}}, \code{\link{additive_penalty}}, \code{\link{distrib_penalty}}
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
#' @description Unlists the aligned theta in the structure's own order.
#' @param pen A \code{StructuredPenalty} object.
#' @param theta The aligned hyperparameter list.
#' @return A numeric vector.
#' @keywords internal
struct_eta <- function(pen, theta) {
  unlist(theta[pen@params], use.names = FALSE)
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
  s <- pen@structure
  eta <- struct_eta(pen, theta)
  om <- parameters7::param_value(s, eta)
  sum(beta * as.numeric(om %*% beta)) / 2 -
    parameters7::param_logdet(s, eta) / 2 + s@rank / 2 * log(2 * pi)
}

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_gradient.StructuredPenalty
#' @keywords internal
S7::method(penalty_gradient, StructuredPenalty) <- function(pen, beta, theta, ...) {
  eta <- struct_eta(pen, theta)
  as.numeric(parameters7::param_value(pen@structure, eta) %*% beta)
}

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_hessian.StructuredPenalty
#' @keywords internal
S7::method(penalty_hessian, StructuredPenalty) <- function(pen, beta, theta, ...) {
  unclass(parameters7::param_value(pen@structure, struct_eta(pen, theta)))
}

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_grad_theta.StructuredPenalty
#' @keywords internal
S7::method(penalty_grad_theta, StructuredPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    s <- pen@structure
    eta <- struct_eta(pen, theta)
    A <- parameters7::param_d1(s, eta)
    dld <- parameters7::param_dlogdet(s, eta)
    stats::setNames(lapply(seq_along(pen@params), function(k) {
      sum(beta * as.numeric(A[[k]] %*% beta)) / 2 - dld[[k]] / 2
    }), pen@params)
  }

#' @rdname penalty_value.StructuredPenalty
#' @name penalty_hess_theta.StructuredPenalty
#' @keywords internal
S7::method(penalty_hess_theta, StructuredPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    s <- pen@structure
    eta <- struct_eta(pen, theta)
    A2 <- parameters7::param_d2(s, eta)
    d2 <- parameters7::param_d2logdet(s, eta)
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
    A <- parameters7::param_d1(pen@structure, struct_eta(pen, theta))
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
  unclass(parameters7::param_value(pen@structure, struct_eta(pen, theta)))
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
  s <- pen@structure
  eta <- struct_eta(pen, theta)
  d1 <- parameters7::param_dlogdet(s, eta)
  d2 <- parameters7::param_d2logdet(s, eta)
  nm <- pen@params
  prs <- ptheta_pairs(nm)
  list(
    value = parameters7::param_logdet(s, eta),
    grad = stats::setNames(as.list(unname(d1)), nm),
    hess = stats::setNames(lapply(prs, function(pr) {
      ij <- sort(match(pr, nm))
      unname(d2[[paste(nm[ij], collapse = ":")]])
    }), names(prs))
  )
}
