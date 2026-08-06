#' @include generics.R
NULL

#' @title S7 Class for the Quadratic Penalty
#'
#' @description
#' The class \code{\link{quadratic_penalty}} instantiates. Beyond the base
#' properties it stores the matrix \eqn{P}, its rank, its null basis and its
#' log pseudo-determinant, all fixed by one eigendecomposition at
#' construction.
#'
#' @inheritParams penalty
#' @param P The symmetric positive semidefinite matrix.
#' @param p_rank The rank of \eqn{P}.
#' @param null_basis An orthonormal basis of the null space of \eqn{D'PD}.
#' @param logpdet_P The log pseudo-determinant of \eqn{P}.
#' @param DPD The assembled \eqn{D'PD}, cached.
#'
#' @return An object of class \code{QuadraticPenalty}.
#'
#' @seealso \code{\link{quadratic_penalty}}
#' @examples
#' S7::S7_inherits(quadratic_penalty(diag(2)), QuadraticPenalty)
#' @keywords internal
#' @export
QuadraticPenalty <- S7::new_class(
  name = "QuadraticPenalty",
  parent = penalty,
  properties = list(
    P = S7::class_any,
    p_rank = S7::class_numeric,
    null_basis = S7::class_any,
    logpdet_P = S7::class_numeric,
    DPD = S7::class_any
  )
)

#' @title Construct a Quadratic Penalty
#'
#' @description
#' The exact negative log-density of the (possibly degenerate) Gaussian
#' prior \eqn{\beta \sim N(0, (\lambda D'PD)^{-})}:
#' \deqn{\rho = \tfrac{\lambda}{2}\,(D\beta)'P(D\beta)
#'   - \tfrac{r}{2}\log\lambda + \tfrac{r}{2}\log 2\pi
#'   - \tfrac{1}{2}\log\mathrm{pdet}(P),}
#' with \eqn{r} the rank of \eqn{P}. The constant is kept deliberately:
#' dropping it, as penalized-likelihood software does, makes joint
#' hyperparameter estimation degenerate.
#'
#' @details
#' One eigendecomposition at construction fixes the rank by the relative
#' rule \eqn{\mathrm{ev} > \mathrm{tol}\cdot\max(\mathrm{ev})} -- a
#' statement about the matrix rather than about whichever arithmetic is
#' later performed on it -- and stores the exact null basis of \eqn{D'PD},
#' so that membership questions never go through a rank recomputed from an
#' assembled sum. Every quantity is then closed form in \eqn{\lambda}.
#'
#' @param P A symmetric positive semidefinite matrix, for instance a
#'   \pkg{basis7} Gram matrix or a difference penalty \eqn{D_k'D_k}.
#' @param map The matrix \eqn{D}, or \code{NULL} (default) for the identity.
#' @param link_lambda The link carrying \eqn{\lambda}; defaults to
#'   \code{linkfunctions7::log_link()}.
#' @param tol The relative eigenvalue tolerance of the rank rule.
#'
#' @return An object of class \code{QuadraticPenalty}.
#'
#' @examples
#' # a second-difference penalty on five coefficients: rank 3, an improper
#' # prior whose null space is the linear polynomials
#' P <- crossprod(diff(diag(5), differences = 2))
#' pen <- quadratic_penalty(P)
#' penalty_rank(pen)
#' is_proper(pen)
#' penalty_value(pen, rnorm(5), list(lambda = 2))
#'
#' @export
quadratic_penalty <- function(P, map = NULL,
                              link_lambda = linkfunctions7::log_link(),
                              tol = 1e-10) {
  P <- as.matrix(P)
  if (nrow(P) != ncol(P) || max(abs(P - t(P))) > 1e-8 * max(1, max(abs(P)))) {
    stop("'P' must be a symmetric matrix.", call. = FALSE)
  }
  P <- (P + t(P)) / 2
  if (!is.null(map)) {
    map <- as.matrix(map)
    if (nrow(map) != nrow(P)) {
      stop("'map' must have as many rows as 'P'.", call. = FALSE)
    }
  }
  q <- if (is.null(map)) nrow(P) else ncol(map)

  ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  keep <- ev > tol * max(ev, 0)
  if (!any(keep)) stop("'P' is the zero matrix.", call. = FALSE)
  r <- sum(keep)
  logpdet_P <- sum(log(ev[keep]))

  DPD <- if (is.null(map)) P else crossprod(map, P %*% map)
  eD <- eigen(DPD, symmetric = TRUE)
  keepD <- eD$values > tol * max(eD$values, 0)
  nb <- eD$vectors[, !keepD, drop = FALSE]

  QuadraticPenalty(
    penalty_name = "quadratic",
    map = map,
    n_coef = q,
    params = "lambda",
    params_bounds = list(lambda = c(0, Inf)),
    link_params = list(lambda = link_lambda),
    params_smooth = c(lambda = TRUE),
    P = P,
    p_rank = r,
    null_basis = nb,
    logpdet_P = logpdet_P,
    DPD = DPD
  )
}

#' The Quadratic Form of a Quadratic Penalty
#' @description The quadratic form of the mapped coefficients, shared by
#'   the value and the
#'   theta derivatives.
#' @param pen A \code{QuadraticPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @return A single number.
#' @keywords internal
quad_form <- function(pen, beta) {
  t <- map_apply(pen, beta)
  sum(t * as.numeric(pen@P %*% t))
}

#' @title Quadratic Penalty Methods
#' @name penalty_value.QuadraticPenalty
#' @description
#' Closed form in \eqn{\lambda} throughout; the Hessian is constant in
#' \eqn{\beta}, and the mixed block is \eqn{D'PD\beta}.
#' @param pen A \code{QuadraticPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A list containing \code{lambda}.
#' @param scale Handled by the generic.
#' @param ... Unused.
#' @return See the generic pages.
#' @keywords internal
S7::method(penalty_value, QuadraticPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta[[1]]
  r <- pen@p_rank
  lam / 2 * quad_form(pen, beta) - r / 2 * log(lam) +
    r / 2 * log(2 * pi) - pen@logpdet_P / 2
}

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_gradient.QuadraticPenalty
#' @keywords internal
S7::method(penalty_gradient, QuadraticPenalty) <- function(pen, beta, theta, ...) {
  theta[[1]] * as.numeric(pen@DPD %*% beta)
}

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_hessian.QuadraticPenalty
#' @keywords internal
S7::method(penalty_hessian, QuadraticPenalty) <- function(pen, beta, theta, ...) {
  theta[[1]] * pen@DPD
}

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_grad_theta.QuadraticPenalty
#' @keywords internal
S7::method(penalty_grad_theta, QuadraticPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    list(lambda = quad_form(pen, beta) / 2 - pen@p_rank / (2 * theta[[1]]))
  }

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_hess_theta.QuadraticPenalty
#' @keywords internal
S7::method(penalty_hess_theta, QuadraticPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    list(lambda_lambda = pen@p_rank / (2 * theta[[1]]^2))
  }

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_cross.QuadraticPenalty
#' @keywords internal
S7::method(penalty_cross, QuadraticPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    list(lambda = as.numeric(pen@DPD %*% beta))
  }

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_kinks.QuadraticPenalty
#' @keywords internal
S7::method(penalty_kinks, QuadraticPenalty) <- function(pen, theta, ...) {
  numeric(0)
}

#' @rdname penalty_value.QuadraticPenalty
#' @name is_proper.QuadraticPenalty
#' @keywords internal
S7::method(is_proper, QuadraticPenalty) <- function(pen, ...) {
  pen@p_rank == ncol(pen@DPD) && ncol(pen@null_basis) == 0L
}

#' @rdname penalty_value.QuadraticPenalty
#' @name is_quadratic.QuadraticPenalty
#' @keywords internal
S7::method(is_quadratic, QuadraticPenalty) <- function(pen, ...) TRUE

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_matrix.QuadraticPenalty
#' @keywords internal
S7::method(penalty_matrix, QuadraticPenalty) <- function(pen, theta, ...) {
  theta[[1]] * pen@DPD
}

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_rank.QuadraticPenalty
#' @keywords internal
S7::method(penalty_rank, QuadraticPenalty) <- function(pen, ...) {
  as.integer(pen@p_rank)
}

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_null_basis.QuadraticPenalty
#' @keywords internal
S7::method(penalty_null_basis, QuadraticPenalty) <- function(pen, ...) {
  pen@null_basis
}

#' @rdname penalty_value.QuadraticPenalty
#' @name penalty_logpdet.QuadraticPenalty
#' @keywords internal
S7::method(penalty_logpdet, QuadraticPenalty) <- function(pen, theta, ...) {
  lam <- theta[[1]]
  r <- pen@p_rank
  list(value = r * log(lam) + pen@logpdet_P,
       grad = list(lambda = r / lam),
       hess = list(lambda_lambda = -r / lam^2))
}
