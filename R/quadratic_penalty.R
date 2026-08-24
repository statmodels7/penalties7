#' @include generics.R
NULL

#' @title S7 Class for the Quadratic Penalty
#'
#' @description
#' The class [quadratic_penalty()] builds. Beyond the seven properties every
#' penalty carries it stores the matrix \eqn{P}, its rank, an orthonormal basis
#' of the null space of \eqn{D'PD}, the log pseudo-determinant of \eqn{P} and
#' the assembled \eqn{D'PD}. All five are fixed by one eigendecomposition at
#' construction and none of them moves with the hyperparameter, so every
#' quantity the branch supplies is closed form in \eqn{\lambda}.
#'
#' @details
#' Construct through [quadratic_penalty()], which computes the five derived
#' properties and checks that \eqn{P} is symmetric. This raw constructor takes
#' them as given and checks nothing, so a wrong rank here would be reported by
#' [penalty_logpdet()] and by nothing else.
#'
#' @inheritParams penalty
#' @param P The symmetric positive semidefinite matrix of the quadratic form,
#'   with as many rows as the map has, or as many as there are coefficients
#'   when the map is `NULL`. A `dgCMatrix` under `blocks > 1`.
#' @param p_rank The rank of \eqn{P}, a single whole number, counted at
#'   construction by a relative eigenvalue rule.
#' @param null_basis An orthonormal basis of the null space of \eqn{D'PD}, with
#'   `n_coef` rows and `n_coef - p_rank` columns, and no columns when the
#'   penalty is full rank.
#' @param logpdet_P The log pseudo-determinant of \eqn{P}: the sum of the
#'   logarithms of its non-zero eigenvalues. A single number.
#' @param DPD The assembled \eqn{D'PD}, an `n_coef` by `n_coef` matrix, cached
#'   so that the gradient and the Hessian never re-form it.
#'
#' @return An S7 object of class `QuadraticPenalty`, inheriting from
#'   [penalty()], with the seven inherited properties and the five above.
#'
#' @seealso [quadratic_penalty()] for the constructor to use,
#'   [penalty_value.QuadraticPenalty()] for what the branch computes,
#'   [structured_penalty()] for the other branch that answers
#'   [is_quadratic()] with `TRUE`.
#'
#' @examples
#' pen <- quadratic_penalty(crossprod(diff(diag(5), differences = 2)))
#' S7::S7_inherits(pen, QuadraticPenalty)
#'
#' # The five derived properties, fixed at construction.
#' pen@p_rank
#' dim(pen@null_basis)
#' pen@logpdet_P
#'
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
#' Builds the penalty \eqn{\tfrac{\lambda}{2}(D\beta)'P(D\beta)} together with
#' the constant that makes it a density, so the value returned is exactly the
#' negative log-density of the Gaussian prior with precision
#' \eqn{\lambda D'PD}. This is the penalty behind a spline smooth, a ridge and
#' a Gaussian random effect: one fixed matrix saying which directions are
#' expensive, and one smoothing parameter saying how expensive.
#'
#' @details
#' # The value
#'
#' With \eqn{r} the rank of \eqn{P},
#'
#' \deqn{\rho(\beta; \lambda) = \tfrac{\lambda}{2}\,(D\beta)'P(D\beta)
#'   - \tfrac{r}{2}\log\lambda + \tfrac{r}{2}\log 2\pi
#'   - \tfrac{1}{2}\log\mathrm{pdet}(P).}
#'
#' The last three terms are the normalizing constant of the prior, taken over
#' the range of \eqn{P} alone when \eqn{P} is deficient. Penalized-likelihood
#' software usually drops them, and dropping them makes \eqn{\lambda}
#' unestimable: with no \eqn{-\tfrac{r}{2}\log\lambda} the penalty falls to
#' zero as \eqn{\lambda} does and the joint maximum runs away.
#'
#' # The rank, and why it is fixed at construction
#'
#' One eigendecomposition fixes the rank by the relative rule
#' \eqn{\mathrm{ev} > \mathtt{tol} \cdot \max(\mathrm{ev})}, and stores the
#' exact null basis of \eqn{D'PD}. That is a statement about the matrix, not
#' about whichever arithmetic is later performed on it, and it is what keeps
#' the answer stable: a rank recounted from an assembled sum of several
#' penalties falls as their smoothing parameters spread apart, while the true
#' null space is the intersection of the components' and does not move.
#'
#' # Blockwise repetition
#'
#' `blocks = m` builds the penalty of \eqn{I_m \otimes P} from \eqn{P} alone,
#' the shape one copy of a smooth per level of a factor takes. The
#' eigenvalues of \eqn{I_m \otimes P} are \eqn{P}'s repeated \eqn{m} times, so
#' the rank is \eqn{mr}, the log pseudo-determinant is
#' \eqn{m\log\mathrm{pdet}(P)} and the null space is \eqn{I_m \otimes N}: the
#' large matrix is never decomposed. Measured on second differences over ten
#' coefficients at \eqn{m = 200}, so an assembled matrix of
#' \eqn{2000 \times 2000}:
#'
#' | | assembled | `blocks = 200` |
#' |---|---|---|
#' | eigendecomposition | 4.79 s | 2.9e-05 s |
#' | whole construction | 5.82 s | 0.001 s |
#' | stored matrix | 30.5 MB dense | 0.11 MB sparse |
#'
#' at a density of 0.0022. The two penalties agree exactly: the value to
#' 1.4e-14, the gradient and the Hessian to 0, and the same rank, null basis
#' and log pseudo-determinant.
#'
#' **Under `blocks > 1` the matrix-valued returns are `dgCMatrix`**, not base
#' matrices: [penalty_hessian()], [penalty_matrix()], [penalty_null_basis()]
#' and [penalty_dhessian()] all keep the sparse storage. A branch built to
#' avoid the assembled form would defeat itself by densifying at its own
#' boundary. Consumers that write these into a block of their own information
#' need to accept both classes.
#'
#' `blocks` does not combine with `map`, and the constructor rejects the pair:
#' a map mixes the blocks, so \eqn{D'(I_m \otimes P)D} is not block diagonal
#' with the structure `blocks` names.
#'
#' @param P A symmetric positive semidefinite matrix, for instance a
#'   \pkg{basis7} Gram matrix or a difference penalty \eqn{D_k'D_k}. Symmetry
#'   is checked to a relative tolerance of `1e-8` and the matrix is symmetrized
#'   before use; a zero matrix is rejected.
#' @param map The matrix \eqn{D}, with as many rows as \eqn{P} and one column
#'   per coefficient, or `NULL` (the default) for the identity, and then
#'   `n_coef` is `nrow(P)`. A \pkg{Matrix} object keeps its own storage.
#' @param blocks How many times \eqn{P} is repeated blockwise, a whole number
#'   of at least 1. `1L`, the default, is the ordinary single-block penalty.
#'   Above 1 it must be given without a `map`.
#' @param link_lambda The \pkg{linkfunctions7} link carrying \eqn{\lambda} onto
#'   the whole real line. `linkfunctions7::log_link()` by default, \eqn{\lambda}
#'   being positive.
#' @param tol The relative eigenvalue tolerance of the rank rule, `1e-10` by
#'   default. An eigenvalue at or below `tol` times the largest counts as zero,
#'   so its direction joins the null space and is not penalized.
#'
#' @return A [QuadraticPenalty()] object with one hyperparameter, `lambda`,
#'   bounded on \eqn{(0, \infty)}.
#'
#' @examples
#' # Second differences over five coefficients: rank 3, and a null space
#' # holding the straight lines, so this is an improper prior.
#' P <- crossprod(diff(diag(5), differences = 2))
#' pen <- quadratic_penalty(P)
#' penalty_rank(pen)
#' is_proper(pen)
#'
#' # A straight line lies in the null space and has zero gradient.
#' penalty_gradient(pen, c(1, 2, 3, 4, 5), list(lambda = 2))
#'
#' # A full-rank penalty is a proper Gaussian prior, and its value is the
#' # negative log-density of one.
#' ridge <- quadratic_penalty(diag(3))
#' b <- c(0.4, -1.1, 0.7)
#' all.equal(penalty_value(ridge, b, list(lambda = 2)),
#'           -sum(stats::dnorm(b, sd = 1 / sqrt(2), log = TRUE)))
#'
#' # A map penalizes what it selects. Here only the second differences of a
#' # six-vector are charged for, through a 4 x 6 map.
#' D <- diff(diag(6), differences = 2)
#' curved <- quadratic_penalty(diag(4), map = D)
#' curved@n_coef
#' penalty_value(curved, 1:6, list(lambda = 1))
#'
#' # blocks = m repeats P without assembling I_m (x) P.
#' blocked <- quadratic_penalty(P, blocks = 4)
#' c(blocked@n_coef, penalty_rank(blocked))
#' all.equal(penalty_rank(blocked), 4L * penalty_rank(pen))
#'
#' @seealso [additive_penalty()] for a sum of quadratics with a smoothing
#'   parameter each, [structured_penalty()] for a matrix that moves with
#'   several hyperparameters, [distrib_penalty()] for the separable branch,
#'   [penalty_logpdet()] for the constant a marginal criterion reads,
#'   [check_penalty()] to verify the result.
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
#' Builds the penalty of \eqn{I_m \otimes P} from \eqn{P} alone, and returns a
#' [QuadraticPenalty()] whose stored matrix is sparse. Called by
#' [quadratic_penalty()] when `blocks` is above one.
#'
#' @details
#' A quadratic penalty needs three things from its matrix: the rank, the log
#' pseudo-determinant and a basis of the null space. All three follow from
#' \eqn{P}, because the eigenvalues of \eqn{I_m \otimes P} are \eqn{P}'s
#' repeated \eqn{m} times: the rank is \eqn{mr}, the log pseudo-determinant is
#' \eqn{m\log\mathrm{pdet}(P)}, and the null space is \eqn{I_m \otimes N}. The
#' large matrix is assembled with [Matrix::kronecker()] and stored, never
#' decomposed.
#'
#' The saving is the whole of the construction. On second differences over ten
#' coefficients at \eqn{m = 200}, the assembled route spends 4.79 s in
#' `eigen()` against 2.9e-05 s here, and 5.82 s in all against 0.001 s. The
#' stored matrix is 0.11 MB against 30.5 MB, at a density of 0.0022, which
#' follows from the sparse storage.
#'
#' The same identity carries [parameters7::kron_identity()] on the other side
#' of the toolkit, for the covariance of grouped random effects.
#'
#' @param P The symmetric matrix of one block, already symmetrized.
#' @param m How many blocks, a whole number of at least two in practice.
#' @param link_lambda The \pkg{linkfunctions7} link for `lambda`.
#' @param tol The relative eigenvalue tolerance for counting a zero.
#'
#' @return A [QuadraticPenalty()] with `n_coef` equal to `m * nrow(P)`, whose
#'   `P`, `DPD` and `null_basis` are `dgCMatrix` objects. Rejects a `P` whose
#'   eigenvalues are all below the tolerance.
#'
#' @seealso [quadratic_penalty()], [parameters7::kron_identity()]
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
#'
#' @description
#' Computes \eqn{(D\beta)'P(D\beta)}, the part of the value that depends on the
#' coefficients. Shared by [penalty_value.QuadraticPenalty()] and the
#' hyperparameter derivatives, both of which need it and neither of which needs
#' anything else from \eqn{\beta}.
#'
#' @details
#' Uses \eqn{P} and the map, not the cached \eqn{D'PD}, so the intermediate is
#' the \eqn{m}-vector \eqn{D\beta} and no \eqn{q \times q} product is formed.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#'
#' @return A single number, non-negative for a positive semidefinite \eqn{P}.
#'
#' @seealso [penalty_value.QuadraticPenalty()]
#'
#' @keywords internal
quad_form <- function(pen, beta) {
  t <- map_apply(pen, beta)
  sum(t * as.numeric(pen@P %*% t))
}

#' @title Value of a Quadratic Penalty
#' @name penalty_value.QuadraticPenalty
#'
#' @description
#' Returns the negative log-density of the Gaussian prior with precision
#' \eqn{\lambda D'PD}, evaluated by adding the constant to the quadratic form.
#' Every term is closed in \eqn{\lambda}, the rank and the log
#' pseudo-determinant of \eqn{P} having been fixed at construction.
#'
#' @details
#' \deqn{\rho(\beta; \lambda) = \frac{\lambda}{2}\,(D\beta)'P(D\beta)
#'   - \frac{r}{2}\log\lambda + \frac{r}{2}\log 2\pi
#'   - \frac{1}{2}\log\mathrm{pdet}(P),}
#'
#' with \eqn{r} the stored rank. At \eqn{\beta = 0} the value is the constant
#' alone and still moves with \eqn{\lambda}, so the smoothing parameter stays
#' estimable from a joint objective.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A list holding `lambda`, already aligned and bound-checked by
#'   the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number.
#'
#' @examples
#' # A full-rank penalty is a proper prior, so the value is a log-density.
#' pen <- quadratic_penalty(diag(3))
#' b <- c(0.4, -1.1, 0.7)
#' all.equal(penalty_value(pen, b, list(lambda = 2)),
#'           -sum(stats::dnorm(b, sd = 1 / sqrt(2), log = TRUE)))
#'
#' # The constant alone, at the origin, and how it moves with lambda.
#' sapply(c(0.5, 2, 8),
#'        function(l) penalty_value(pen, c(0, 0, 0), list(lambda = l)))
#'
#' @seealso [penalty_gradient.QuadraticPenalty()] for the coefficient
#'   derivatives, [penalty_logpdet.QuadraticPenalty()] for the constant's own
#'   piece, [quadratic_penalty()] for the constructor.
#' @keywords internal
S7::method(penalty_value, QuadraticPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta[[1]]
  r <- pen@p_rank
  lam / 2 * quad_form(pen, beta) - r / 2 * log(lam) +
    r / 2 * log(2 * pi) - pen@logpdet_P / 2
}

#' @title Coefficient Derivatives of a Quadratic Penalty
#' @name penalty_gradient.QuadraticPenalty
#'
#' @description
#' `penalty_gradient()` returns \eqn{\lambda D'PD\beta} and `penalty_hessian()`
#' returns \eqn{\lambda D'PD}. Both read the cached \eqn{D'PD}, so neither
#' touches the map at call time, and the Hessian does not depend on \eqn{\beta}
#' at all, the value being a quadratic form.
#'
#' @details
#' Because the Hessian is constant, the gradient is exactly the Hessian applied
#' to the coefficients, and a Newton step on this penalty alone reaches the
#' minimum in one iteration from anywhere.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A list holding `lambda`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` a numeric vector of length `pen@n_coef`;
#'   `penalty_hessian()` a symmetric matrix of that side, a base matrix for an
#'   ordinary penalty and a `dgCMatrix` when the penalty was built with
#'   `blocks > 1`.
#'
#' @examples
#' pen <- quadratic_penalty(crossprod(diff(diag(4))))
#' b <- c(1, -0.5, 0.2, 2)
#'
#' # The gradient is the constant Hessian applied to the coefficients.
#' all.equal(penalty_gradient(pen, b, list(lambda = 3)),
#'           drop(penalty_hessian(pen, b, list(lambda = 3)) %*% b))
#'
#' # The Hessian does not move with beta.
#' all.equal(penalty_hessian(pen, b, list(lambda = 3)),
#'           penalty_hessian(pen, 0 * b, list(lambda = 3)))
#'
#' @seealso [penalty_value.QuadraticPenalty()] for the quantity
#'   differentiated, [penalty_grad_theta.QuadraticPenalty()] for the
#'   hyperparameter derivatives.
#' @keywords internal
S7::method(penalty_gradient, QuadraticPenalty) <- function(pen, beta, theta, ...) {
  theta[[1]] * as.numeric(pen@DPD %*% beta)
}

#' @rdname penalty_gradient.QuadraticPenalty
#' @name penalty_hessian.QuadraticPenalty
#' @keywords internal
S7::method(penalty_hessian, QuadraticPenalty) <- function(pen, beta, theta, ...) {
  theta[[1]] * pen@DPD
}

#' @title Hyperparameter Derivatives of a Quadratic Penalty
#' @name penalty_grad_theta.QuadraticPenalty
#'
#' @description
#' The three blocks in \eqn{\lambda}, all closed form. `penalty_grad_theta()`
#' returns \eqn{\tfrac{1}{2}(D\beta)'P(D\beta) - r/(2\lambda)},
#' `penalty_hess_theta()` returns \eqn{r/(2\lambda^2)}, and `penalty_cross()`
#' returns \eqn{D'PD\beta}, which is the gradient with \eqn{\lambda} divided
#' out.
#'
#' @details
#' The value is affine in \eqn{\lambda} apart from the
#' \eqn{-\tfrac{r}{2}\log\lambda} of the normalizing constant, so the second
#' derivative comes from that term alone and is positive: the value is convex
#' in the smoothing parameter, with its minimum at
#' \eqn{\lambda = r / (D\beta)'P(D\beta)}. That is the joint-mode estimate of
#' \eqn{\lambda} at fixed coefficients, and it is finite only because the
#' constant was kept.
#'
#' The mixed block carries no \eqn{\lambda} at all, the value being linear in
#' it in every coefficient direction.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A list holding `lambda`.
#' @param scale Read by the generic, which applies the chain rule onto the
#'   unconstrained scale after dispatch. These methods always return the
#'   parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` a list of one number named `lambda`;
#'   `penalty_hess_theta()` a list of one number named `lambda_lambda`;
#'   `penalty_cross()` a list of one numeric vector of length `pen@n_coef`,
#'   named `lambda`.
#'
#' @examples
#' pen <- quadratic_penalty(diag(3))
#' b <- c(1, -0.5, 0.2)
#'
#' # The value is convex in lambda, with its minimum where the gradient
#' # vanishes: at r divided by the quadratic form.
#' lam_hat <- 3 / sum(b^2)
#' penalty_grad_theta(pen, b, list(lambda = lam_hat))
#' penalty_hess_theta(pen, b, list(lambda = lam_hat))
#'
#' # The mixed block is the gradient with lambda divided out.
#' all.equal(penalty_cross(pen, b, list(lambda = 4))$lambda,
#'           penalty_gradient(pen, b, list(lambda = 4)) / 4)
#'
#' @seealso [penalty_value.QuadraticPenalty()] for the quantity
#'   differentiated, [penalty_dhessian.QuadraticPenalty()] for the third-order
#'   blocks, [penalty_logpdet.QuadraticPenalty()] for the constant's
#'   derivatives.
#' @keywords internal
S7::method(penalty_grad_theta, QuadraticPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    list(lambda = quad_form(pen, beta) / 2 - pen@p_rank / (2 * theta[[1]]))
  }

#' @rdname penalty_grad_theta.QuadraticPenalty
#' @name penalty_hess_theta.QuadraticPenalty
#' @keywords internal
S7::method(penalty_hess_theta, QuadraticPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    list(lambda_lambda = pen@p_rank / (2 * theta[[1]]^2))
  }

#' @rdname penalty_grad_theta.QuadraticPenalty
#' @name penalty_cross.QuadraticPenalty
#' @keywords internal
S7::method(penalty_cross, QuadraticPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    list(lambda = as.numeric(pen@DPD %*% beta))
  }

#' @title Smoothness and Kind of a Quadratic Penalty
#' @name penalty_kinks.QuadraticPenalty
#'
#' @description
#' The three questions a consumer asks before choosing how to fit a quadratic
#' block. `penalty_kinks()` returns `numeric(0)`, the value being a polynomial
#' and smooth everywhere. `is_quadratic()` returns `TRUE`, so the four marginal
#' quantities are available. `is_proper()` returns `TRUE` only when the penalty
#' has full rank.
#'
#' @details
#' Properness is a rank test and nothing else: \eqn{\exp(-\rho)} is flat along
#' every null direction of \eqn{D'PD}, so it integrates only when there are
#' none. A ridge over \eqn{q} coefficients is proper; second differences over
#' \eqn{q} leave the level and the slope free, have rank \eqn{q - 2}, and are
#' not. Both are usable, and the difference shows in the normalizing constant,
#' which is taken over the range alone.
#'
#' Because there are no kinks, a quadratic block goes to whatever smooth method
#' the fit uses, and [penalty_prox()] is a linear solve rather than a
#' threshold.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param theta A list holding `lambda`. Read by `penalty_kinks()`, whose
#'   answer does not depend on it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` a numeric vector of length zero.
#'   `is_quadratic()` the single logical `TRUE`.
#'   `is_proper()` a single logical, `TRUE` when the rank equals the number of
#'   coefficients and the null basis has no columns.
#'
#' @examples
#' ridge <- quadratic_penalty(diag(4))
#' curve <- quadratic_penalty(crossprod(diff(diag(4), differences = 2)))
#'
#' penalty_kinks(ridge, list(lambda = 1))
#' c(is_quadratic(ridge), is_quadratic(curve))
#' c(is_proper(ridge), is_proper(curve))
#' c(penalty_rank(ridge), penalty_rank(curve))
#'
#' @seealso [penalty_kinks()], [is_proper()] and [is_quadratic()] for the
#'   generics, [penalty_matrix.QuadraticPenalty()] for the quantities
#'   `is_quadratic()` gates.
#' @keywords internal
S7::method(penalty_kinks, QuadraticPenalty) <- function(pen, theta, ...) {
  numeric(0)
}

#' @rdname penalty_kinks.QuadraticPenalty
#' @name is_proper.QuadraticPenalty
#' @keywords internal
S7::method(is_proper, QuadraticPenalty) <- function(pen, ...) {
  pen@p_rank == ncol(pen@DPD) && ncol(pen@null_basis) == 0L
}

#' @rdname penalty_kinks.QuadraticPenalty
#' @name is_quadratic.QuadraticPenalty
#' @keywords internal
S7::method(is_quadratic, QuadraticPenalty) <- function(pen, ...) TRUE

#' @title Marginal Quantities of a Quadratic Penalty
#' @name penalty_matrix.QuadraticPenalty
#'
#' @description
#' The four pieces a REML or marginal likelihood criterion reads.
#' `penalty_matrix()` returns \eqn{\lambda D'PD}, `penalty_rank()` the rank
#' fixed at construction, `penalty_null_basis()` the stored orthonormal basis
#' of the null space, and `penalty_logpdet()` the log pseudo-determinant with
#' its first two derivatives in \eqn{\lambda}.
#'
#' @details
#' Since \eqn{S(\lambda) = \lambda D'PD} and only the range contributes,
#'
#' \deqn{\log^{+}\lvert S(\lambda) \rvert = r\log\lambda
#'   + \log\mathrm{pdet}(P), \qquad
#'   \frac{\partial}{\partial\lambda} = \frac{r}{\lambda}, \qquad
#'   \frac{\partial^2}{\partial\lambda^2} = -\frac{r}{\lambda^2},}
#'
#' all three exact and none of them needing a decomposition at call time. The
#' rank and the null basis do not depend on \eqn{\lambda} and their generics
#' take no `theta`.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param theta A list holding `lambda`. `penalty_matrix()` and
#'   `penalty_logpdet()` only.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_matrix()` a symmetric matrix of side `pen@n_coef`, a base
#'   matrix ordinarily and a `dgCMatrix` when the penalty was built with
#'   `blocks > 1`.
#'   `penalty_rank()` a single integer.
#'   `penalty_null_basis()` a matrix with `pen@n_coef` rows and
#'   `pen@n_coef - penalty_rank(pen)` orthonormal columns, sparse under
#'   `blocks > 1`.
#'   `penalty_logpdet()` a list of `value` (a number), `grad` (a list of one
#'   number, `lambda`) and `hess` (a list of one number, `lambda_lambda`).
#'
#' @examples
#' pen <- quadratic_penalty(crossprod(diff(diag(5), differences = 2)))
#' penalty_rank(pen)
#'
#' # The null basis spans the straight lines, and the matrix kills it.
#' N <- penalty_null_basis(pen)
#' dim(N)
#' max(abs(penalty_matrix(pen, list(lambda = 3)) %*% N))
#'
#' # The log pseudo-determinant is r log(lambda) plus a constant.
#' lp <- penalty_logpdet(pen, list(lambda = 3))
#' c(lp$grad$lambda, lp$hess$lambda_lambda)
#' c(penalty_rank(pen) / 3, -penalty_rank(pen) / 9)
#'
#' @seealso [penalty_matrix()] and its three siblings for the generics,
#'   [is_quadratic()] for the predicate that gates them,
#'   [penalty_value.QuadraticPenalty()] for the value that carries the same
#'   constant.
#' @keywords internal
S7::method(penalty_matrix, QuadraticPenalty) <- function(pen, theta, ...) {
  theta[[1]] * pen@DPD
}

#' @rdname penalty_matrix.QuadraticPenalty
#' @name penalty_rank.QuadraticPenalty
#' @keywords internal
S7::method(penalty_rank, QuadraticPenalty) <- function(pen, ...) {
  as.integer(pen@p_rank)
}

#' @rdname penalty_matrix.QuadraticPenalty
#' @name penalty_null_basis.QuadraticPenalty
#' @keywords internal
S7::method(penalty_null_basis, QuadraticPenalty) <- function(pen, ...) {
  pen@null_basis
}

#' @rdname penalty_matrix.QuadraticPenalty
#' @name penalty_logpdet.QuadraticPenalty
#' @keywords internal
S7::method(penalty_logpdet, QuadraticPenalty) <- function(pen, theta, ...) {
  lam <- theta[[1]]
  r <- pen@p_rank
  list(value = r * log(lam) + pen@logpdet_P,
       grad = list(lambda = r / lam),
       hess = list(lambda_lambda = -r / lam^2))
}
