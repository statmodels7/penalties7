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
#' @param blocks How many times \eqn{P} is repeated blockwise: the penalty is
#'   then that of \eqn{I_m \otimes P}, which is what one copy of a smooth per
#'   level of a factor needs. The big matrix is NEVER formed or decomposed --
#'   the eigenvalues of \eqn{I_m \otimes P} are \eqn{P}'s repeated \eqn{m}
#'   times, so the rank is \eqn{m\,r}, the log pseudo-determinant is
#'   \eqn{m\log\mathrm{pdet}(P)} and the null space is \eqn{I_m \otimes N}.
#'   Measured at \eqn{m = 200} over a basis of ten, that is 4.50 seconds of
#'   eigendecomposition saved and a stored matrix of 0.13 MB against 25.9. It
#'   does not combine with \code{map}, which would mix the blocks.
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
#' @seealso \code{\link{additive_penalty}}, \code{\link{distrib_penalty}}, \code{\link{structured_penalty}}
#' @export
quadratic_penalty <- function(P, map = NULL, blocks = 1L,
                              link_lambda = linkfunctions7::log_link(),
                              tol = 1e-10) {
  P <- as.matrix(P)
  if (nrow(P) != ncol(P) || max(abs(P - t(P))) > 1e-8 * max(1, max(abs(P)))) {
    stop("'P' must be a symmetric matrix.", call. = FALSE)
  }
  P <- (P + t(P)) / 2
  if (!is.numeric(blocks) || length(blocks) != 1L || is.na(blocks) ||
      blocks < 1 || blocks != round(blocks)) {
    stop("'blocks' must be a whole number of at least 1.", call. = FALSE)
  }
  blocks <- as.integer(blocks)
  if (blocks > 1L) {
    if (!is.null(map)) {
      stop(paste0("'blocks' and 'map' do not combine: a map mixes the",
                  " blocks, and\n  D'(I (x) P)D is block diagonal with a",
                  " DIFFERENT block each, which is not\n  the structure",
                  " 'blocks' names."), call. = FALSE)
    }
    return(.kron_quadratic(P, blocks, link_lambda, tol))
  }
  if (!is.null(map)) {
  # A map that is already a Matrix is KEPT as it is: `as.matrix()` here would
  # densify a diagonal or sparse map, which is the whole cost the map exists
  # to avoid -- a diagonal one is a per-coordinate rescaling and costs q
  # numbers, its dense form q^2.
    map <- as_map(map)
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

  # the congruence carries a Matrix map's class into the stored matrix, and
  # the stored matrix is dense at any width in the identity branch already;
  # kept a base matrix so a consumer meets one contract
  DPD <- if (is.null(map)) P else as.matrix(crossprod(map, P %*% map))
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

#' A Quadratic Penalty Repeated Blockwise, Without Forming It
#'
#' @description
#' The penalty of \eqn{I_m \otimes P}, built from \eqn{P} alone.
#'
#' @details
#' What the constructor needs from the matrix is its rank, its log pseudo-
#' determinant and a basis of its null space, and all three of them follow
#' from \eqn{P}: the eigenvalues of \eqn{I_m \otimes P} are \eqn{P}'s
#' repeated \eqn{m} times, so the rank is \eqn{m\,r}, the log pseudo-
#' determinant is \eqn{m\log\mathrm{pdet}(P)} and the null space is
#' \eqn{I_m \otimes N}. The big matrix is therefore never decomposed.
#'
#' It is what one smooth per level of a factor needs, and the saving is the
#' whole of the construction: at \eqn{m = 200} over a basis of ten, the
#' eigendecomposition of the assembled \eqn{1800 \times 1800} matrix costs
#' 4.50 seconds and the one of \eqn{P} costs nothing measurable. The stored
#' matrix is sparse besides -- 25.9 MB dense at a density of 0.0005 -- which
#' follows rather than being the point.
#'
#' The same identity is what \code{\link[parameters7]{kron_identity}} uses on
#' the other side of the toolkit, for the covariance of grouped random
#' effects.
#'
#' @param P The symmetric matrix of one block.
#' @param m How many blocks.
#' @param link_lambda The link for the hyperparameter.
#' @param tol The relative tolerance for a zero eigenvalue.
#'
#' @return A \code{QuadraticPenalty}.
#'
#' @seealso \code{\link{quadratic_penalty}}
#'
#' @keywords internal
.kron_quadratic <- function(P, m, link_lambda, tol) {
  e <- eigen(P, symmetric = TRUE)
  keep <- e$values > tol * max(e$values, 0)
  if (!any(keep)) stop("'P' is the zero matrix.", call. = FALSE)
  I <- Matrix::Diagonal(m)
  Pb <- Matrix::kronecker(I, methods::as(Matrix::Matrix(P, sparse = TRUE),
                                         "generalMatrix"))
  nb <- Matrix::kronecker(I, Matrix::Matrix(e$vectors[, !keep, drop = FALSE],
                                            sparse = TRUE))
  QuadraticPenalty(
    penalty_name = "quadratic",
    map = NULL,
    n_coef = m * nrow(P),
    params = "lambda",
    params_bounds = list(lambda = c(0, Inf)),
    link_params = list(lambda = link_lambda),
    params_smooth = c(lambda = TRUE),
    P = Pb,
    p_rank = m * sum(keep),
    null_basis = nb,
    logpdet_P = m * sum(log(e$values[keep])),
    DPD = Pb
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
