#' @include generics.R
NULL

#' @title S7 Class for the Quadratic Penalty
#'
#' @description
#' `QuadraticPenalty` is the class of the objects that [quadratic_penalty()]
#' builds. Beyond the seven properties that every penalty carries, it stores
#' the matrix \eqn{P}, the rank of \eqn{D'PD}, an orthonormal basis of the null
#' space of \eqn{D'PD}, the log pseudo-determinant of \eqn{D'PD} and the
#' assembled \eqn{D'PD}. All five are computed once at construction and are
#' independent of the hyperparameter, so every quantity of the branch is closed
#' form in \eqn{\lambda}.
#'
#' @details
#' Construct through [quadratic_penalty()], which computes the five
#' properties and checks that \eqn{P} is symmetric and positive semidefinite.
#' This raw constructor takes them as given and does not check any of them. A
#' wrong rank passed here enters the value, the hyperparameter derivatives,
#' [penalty_rank()], [penalty_logpdet()] and [is_proper()], and
#' [check_penalty()] does not detect it, because every one of its rows reads
#' the same stored rank.
#'
#' @section Methods:
#' This section lists the methods registered on this class, by the page that
#' documents each.
#'
#' * The value and its coefficient derivatives:
#'   [penalty_value()][penalty_value.QuadraticPenalty],
#'   [penalty_gradient()][penalty_gradient.QuadraticPenalty] and
#'   [penalty_hessian()][penalty_hessian.QuadraticPenalty].
#' * The hyperparameter derivatives:
#'   [penalty_grad_theta()][penalty_grad_theta.QuadraticPenalty],
#'   [penalty_hess_theta()][penalty_hess_theta.QuadraticPenalty] and
#'   [penalty_cross()][penalty_cross.QuadraticPenalty].
#' * The predicates: [penalty_kinks()][penalty_kinks.QuadraticPenalty],
#'   [is_proper()][is_proper.QuadraticPenalty] and
#'   [is_quadratic()][is_quadratic.QuadraticPenalty].
#' * The quantities of a marginal criterion:
#'   [penalty_matrix()][penalty_matrix.QuadraticPenalty],
#'   [penalty_rank()][penalty_rank.QuadraticPenalty],
#'   [penalty_null_basis()][penalty_null_basis.QuadraticPenalty] and
#'   [penalty_logpdet()][penalty_logpdet.QuadraticPenalty].
#' * The higher derivatives of a marginal criterion:
#'   [penalty_dhessian()][penalty_dhessian.QuadraticPenalty],
#'   [penalty_d2hessian()][penalty_d2hessian.QuadraticPenalty],
#'   [penalty_dcross()][penalty_dcross.QuadraticPenalty] and
#'   [beta_quadratic()][beta_quadratic.QuadraticPenalty].
#' * The movement of the Hessian with the coefficients:
#'   [penalty_dhessian_beta()][penalty_dhessian_beta.QuadraticPenalty],
#'   [penalty_d2hessian_beta()][penalty_d2hessian_beta.penalty] and
#'   [penalty_dhessian_beta_theta()][penalty_dhessian_beta_theta.penalty].
#'   The last two return zero matrices because the Hessian does not depend on
#'   the coefficients, and are documented on the page of the base-class
#'   methods.
#' * The proximal operator:
#'   [penalty_prox()][penalty_prox.QuadraticPenalty].
#'
#' The class has no [penalty_prox_spec()], [penalty_draw()] or
#' [penalty_readable()] method of its own, so those generics dispatch to the
#' methods of the base class [penalty()].
#'
#' @inheritParams penalty
#' @param P The symmetric positive semidefinite matrix of the quadratic form,
#'   with as many rows as the map has, or as many as there are coefficients
#'   when the map is `NULL`. It is a `dgCMatrix` under `blocks > 1`.
#' @param p_rank The rank of \eqn{D'PD}, a single whole number, counted at
#'   construction by a relative eigenvalue rule. Without a map it is the rank
#'   of \eqn{P}.
#' @param null_basis An orthonormal basis of the null space of \eqn{D'PD}, with
#'   `n_coef` rows and `n_coef - p_rank` columns, and no columns when the
#'   penalty is full rank.
#' @param logpdet_DPD The log pseudo-determinant of \eqn{D'PD}: the sum of the
#'   logarithms of its non-zero eigenvalues. It is a single number. Without a
#'   map it is the log pseudo-determinant of \eqn{P}.
#' @param DPD The assembled \eqn{D'PD}, an `n_coef` by `n_coef` matrix, cached
#'   so that the gradient and the Hessian never re-form it.
#' @param tol The relative eigenvalue tolerance, a single positive number,
#'   `1e-10` by default. [quadratic_penalty()] counts `p_rank` at this
#'   tolerance and stores it, and [check_penalty()] counts the eigenvalues of
#'   the matrix at the same tolerance when it compares them with the rank.
#'
#' @return An S7 object of class `QuadraticPenalty`, inheriting from
#'   [penalty()], with the seven inherited properties and the six above.
#'
#' @seealso [quadratic_penalty()] for the constructor to use,
#'   [penalty_value.QuadraticPenalty()] for what the branch computes,
#'   [structured_penalty()] and [additive_penalty()] for the other branches
#'   for which [is_quadratic()] returns `TRUE`.
#'
#' @examples
#' pen <- quadratic_penalty(crossprod(diff(diag(5), differences = 2)))
#' S7::S7_inherits(pen, QuadraticPenalty)
#'
#' # These are three of the six properties, all computed at construction.
#' pen@p_rank
#' dim(pen@null_basis)
#' pen@logpdet_DPD
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
    logpdet_DPD = S7::class_numeric,
    DPD = S7::class_any,
    tol = S7::new_property(S7::class_numeric, default = 1e-10)
  )
)

#' @title Construct a Quadratic Penalty
#'
#' @description
#' Builds the penalty \eqn{\tfrac{\lambda}{2}(D\beta)'P(D\beta)} together with
#' the constant that makes it a density, so the value returned is exactly the
#' negative log-density of the Gaussian prior with precision
#' \eqn{\lambda D'PD}. This is the penalty of a spline smooth, a ridge and a
#' Gaussian random effect: one fixed matrix that weights the directions of the
#' coefficients, and one smoothing parameter that scales the weights.
#'
#' @details
#' # The value
#'
#' With \eqn{r} the rank of \eqn{D'PD},
#'
#' \deqn{\rho(\beta; \lambda) = \tfrac{\lambda}{2}\,(D\beta)'P(D\beta)
#'   - \tfrac{r}{2}\log\lambda + \tfrac{r}{2}\log 2\pi
#'   - \tfrac{1}{2}\log\mathrm{pdet}(D'PD).}
#'
#' The last three terms are the normalizing constant of the prior, taken only
#' over the range of \eqn{D'PD} when that matrix is rank deficient. The
#' constant is the one of the matrix applied to \eqn{\beta}, so it carries the
#' map. For a diagonal map with non-zero entries \eqn{d_j} and a full-rank
#' \eqn{P},
#' \eqn{\log\mathrm{pdet}(D'PD) = \log\det P + 2\sum_j\log\lvert d_j\rvert},
#' and \eqn{\sum_j\log\lvert d_j\rvert} is the log Jacobian of
#' \eqn{\beta \mapsto D\beta}. The value is then the negative log-density of a
#' proper Gaussian prior on \eqn{\beta}, and a marginal criterion built on it
#' does not change when a covariate is multiplied by a factor and the map
#' multiplies the coefficient of that covariate by the same factor. The
#' constant also makes \eqn{\lambda} estimable from a joint objective: without
#' the term
#' \eqn{-\tfrac{r}{2}\log\lambda} the penalty decreases to zero with
#' \eqn{\lambda}, and the joint maximum moves to the boundary
#' \eqn{\lambda \to 0}.
#'
#' # The rank and the null basis
#'
#' The rank is counted at construction by the relative rule
#' \eqn{\mathrm{ev} > \mathtt{tol} \cdot \max(\mathrm{ev})} on the eigenvalues
#' of \eqn{D'PD}, and the same mask selects the orthonormal null basis among
#' the eigenvectors, so the null basis always has \eqn{q - r} columns, with
#' \eqn{q} the number of coefficients. Without
#' a map the eigenvalues come from a values-only decomposition of \eqn{P} and
#' the eigenvectors from a second decomposition of the same matrix; both are
#' sorted in decreasing order, so the mask selects the matching columns. With
#' a map the rank, the log pseudo-determinant and the null basis all come from
#' the decomposition of \eqn{D'PD}.
#'
#' Fixing the rank on each matrix keeps it stable when quadratic penalties are
#' combined. A rank recounted from an assembled sum of several penalties falls
#' as their smoothing parameters spread apart, while the true null space is the
#' intersection of the null spaces of the components and does not depend on
#' the smoothing parameters.
#'
#' The rule is relative to the largest eigenvalue of \eqn{D'PD} itself. A map
#' that annihilates \eqn{P} only up to rounding, so that every entry of
#' \eqn{D'PD} is a rounding error relative to \eqn{P}, is therefore not
#' detected: the largest of those rounding errors passes the rule, and the
#' penalty is built with a positive rank and a log pseudo-determinant made of
#' rounding errors. The constructor rejects the map only when \eqn{D'PD} has no
#' positive eigenvalue, as when \eqn{D'PD} is exactly the zero matrix.
#'
#' # Blockwise repetition
#'
#' `blocks = m` builds the penalty of \eqn{I_m \otimes P} from \eqn{P} alone,
#' which is the form that one copy of a smooth per level of a factor takes. The
#' eigenvalues of \eqn{I_m \otimes P} are those of \eqn{P} repeated \eqn{m}
#' times, so the rank is \eqn{mr}, the log pseudo-determinant is
#' \eqn{m\log\mathrm{pdet}(P)} and the null space is \eqn{I_m \otimes N}, with
#' \eqn{N} the null basis of \eqn{P}. Only \eqn{P} is decomposed, and the
#' large matrix is stored in sparse form. The result has the same value,
#' gradient, Hessian, rank, null basis and log pseudo-determinant as the
#' penalty built from the assembled \eqn{I_m \otimes P}, up to rounding.
#'
#' Under `blocks > 1` the matrix-valued returns are `dgCMatrix`, not base
#' matrices: [penalty_hessian()], [penalty_matrix()], [penalty_null_basis()]
#' and the matrices of [penalty_dhessian()] all keep the sparse storage, so the
#' assembled matrix is never densified. A consumer that writes these results
#' into a block of its own information matrix has to accept both classes.
#'
#' `blocks` does not combine with `map`, and the constructor rejects the pair:
#' a map mixes the blocks, so \eqn{D'(I_m \otimes P)D} is in general not block
#' diagonal and does not have the structure that `blocks` builds.
#'
#' @param P A symmetric positive semidefinite matrix, for instance a
#'   \pkg{basis7} Gram matrix or a difference penalty \eqn{D_k'D_k}. A matrix
#'   with an entry that is not finite (`NA`, `NaN` or infinite) is rejected.
#'   Symmetry is checked to a tolerance of `1e-8` times the larger of 1 and
#'   the largest absolute entry of `P`, and the matrix is symmetrized before
#'   use. A matrix whose largest eigenvalue is not positive (a zero or
#'   negative semidefinite matrix) is rejected, and so is one with an
#'   eigenvalue below `-tol` times the largest.
#' @param map The matrix \eqn{D}, with as many rows as \eqn{P} and one column
#'   per coefficient, or `NULL` (the default) for the identity, and then
#'   `n_coef` is `nrow(P)`. A \pkg{Matrix} object keeps its own storage. A map
#'   with an entry that is not finite is rejected, and so is a map for which
#'   \eqn{D'PD} has no positive eigenvalue; a map that
#'   annihilates \eqn{P} only up to rounding is not detected (see the
#'   details).
#' @param blocks How many times \eqn{P} is repeated blockwise, a whole number
#'   of at least 1. `1L`, the default, is the ordinary single-block penalty.
#'   Above 1 it must be given without a `map`.
#' @param link_lambda The \pkg{linkfunctions7} link carrying \eqn{\lambda} onto
#'   the whole real line. The default is `linkfunctions7::log_link()`, because
#'   \eqn{\lambda} is positive.
#' @param tol The relative eigenvalue tolerance of the rank rule, `1e-10` by
#'   default. An eigenvalue of \eqn{D'PD} at or below `tol` times the largest
#'   counts as zero, so its direction joins the null space and is not
#'   penalized. The same tolerance bounds the negative eigenvalues of \eqn{P}
#'   that are accepted as rounding errors. The object stores `tol`, and
#'   [check_penalty()] counts the rank at the same tolerance.
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
#' # A map penalizes its image. Here only the second differences of a
#' # six-vector are penalized, through a 4 x 6 map.
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
#'   [penalty_logpdet()] for the constant that a marginal criterion reads,
#'   [check_penalty()] to verify the result.
#' @export
quadratic_penalty <- function(P, map = NULL, blocks = 1L,
                              link_lambda = linkfunctions7::log_link(),
                              tol = 1e-10) {
  P <- as.matrix(P)
  if (!all(is.finite(P))) {
    stop("'P' must have only finite entries.", call. = FALSE)
  }
  if (nrow(P) != ncol(P) || max(abs(P - t(P))) > 1e-8 * max(1, max(abs(P)))) {
    stop("'P' must be a symmetric matrix.", call. = FALSE)
  }
  P <- (P + t(P)) / 2
  if (!is.numeric(blocks) || length(blocks) != 1L || is.na(blocks) ||
      blocks < 1 || blocks != round(blocks)) {
    stop("'blocks' must be a whole number of at least 1.", call. = FALSE)
  }
  blocks <- as.integer(blocks)
  # the same two checks additive_penalty() applies to each component
  ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
  if (max(ev) <= 0) {
    stop("'P' is zero or negative semidefinite.", call. = FALSE)
  }
  if (min(ev) < -tol * max(ev)) {
    stop("'P' is not positive semidefinite.", call. = FALSE)
  }
  if (blocks > 1L) {
    if (!is.null(map)) {
      stop(paste0("'blocks' and 'map' do not combine: a map mixes the",
                  " blocks, so D'(I (x) P)D\n  is in general not block",
                  " diagonal and does not have the structure that\n",
                  "  'blocks' builds."), call. = FALSE)
    }
    return(.kron_quadratic(P, blocks, link_lambda, tol))
  }
  if (!is.null(map)) {
  # A map that is already a Matrix is KEPT as it is: `as.matrix()` here would
  # densify a diagonal or sparse map, which is the whole cost the map exists
  # to avoid -- a diagonal one is a per-coordinate rescaling and costs q
  # numbers, its dense form q^2.
    map <- as_map(map)
    # a Matrix map is read through its stored entries, never densified; one
    # without an `x` slot is a pattern matrix, whose entries are ones
    if (!all(is.finite(if (!isS4(map)) map
                       else if (methods::.hasSlot(map, "x")) map@x))) {
      stop("'map' must have only finite entries.", call. = FALSE)
    }
    if (nrow(map) != nrow(P)) {
      stop("'map' must have as many rows as 'P'.", call. = FALSE)
    }
  }
  q <- if (is.null(map)) nrow(P) else ncol(map)

  keep <- ev > tol * max(ev, 0)
  r <- sum(keep)
  logpdet <- sum(log(ev[keep]))

  # the congruence carries a Matrix map's class into the stored matrix, and
  # the stored matrix is dense at any width in the identity branch already;
  # kept a base matrix so a consumer meets one contract
  DPD <- if (is.null(map)) P else as.matrix(crossprod(map, P %*% map))
  eD <- eigen(DPD, symmetric = TRUE)
  # With a map the constant is the one of D'PD, the matrix applied to beta, so
  # that it carries the Jacobian of beta -> D beta. Without one the two
  # matrices are the same, and the values-only decomposition above is kept so
  # that every penalty without a map is unchanged to the last bit; its mask
  # also selects the null basis, both decompositions being sorted decreasing,
  # so that the rank and the null basis are read from the same eigenvalues.
  if (is.null(map)) {
    nb <- eD$vectors[, !keep, drop = FALSE]
  } else {
    keepD <- eD$values > tol * max(eD$values, 0)
    if (!any(keepD)) {
      stop("'map' annihilates 'P': D'PD is the zero matrix.", call. = FALSE)
    }
    r <- sum(keepD)
    logpdet <- sum(log(eD$values[keepD]))
    nb <- eD$vectors[, !keepD, drop = FALSE]
  }

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
    logpdet_DPD = logpdet,
    DPD = DPD,
    tol = tol
  )
}

#' A Quadratic Penalty Repeated Blockwise, Without Forming It
#'
#' @description
#' Builds the penalty of \eqn{I_m \otimes P} from \eqn{P} alone, and returns a
#' [QuadraticPenalty()] whose stored matrix is sparse.
#' [quadratic_penalty()] calls it when `blocks` is above one.
#'
#' @details
#' Three quantities of a quadratic penalty come from its matrix: the rank, the
#' log pseudo-determinant and a basis of the null space. All three follow from
#' \eqn{P}, because the eigenvalues of \eqn{I_m \otimes P} are those of
#' \eqn{P} repeated \eqn{m} times: the rank is \eqn{mr}, the log
#' pseudo-determinant is \eqn{m\log\mathrm{pdet}(P)}, and the null space is
#' \eqn{I_m \otimes N}. One eigendecomposition of \eqn{P} supplies the
#' eigenvalues and the eigenvectors, and one mask selects both the rank and
#' the null basis. The large matrix is assembled in sparse form with
#' [Matrix::kronecker()] and stored, and it is never decomposed, so the cost
#' of the construction is that of decomposing \eqn{P}.
#'
#' The same identity underlies [parameters7::kron_identity()], for the
#' covariance of grouped random effects.
#'
#' @param P The symmetric matrix of one block, already symmetrized.
#' @param m How many blocks, a whole number of at least two in practice.
#' @param link_lambda The \pkg{linkfunctions7} link for `lambda`.
#' @param tol The relative eigenvalue tolerance for counting a zero.
#'
#' @return A [QuadraticPenalty()] with `n_coef` equal to `m * nrow(P)`, whose
#'   `P`, `DPD` and `null_basis` are `dgCMatrix` objects. `P` is checked to be
#'   positive semidefinite and not zero by [quadratic_penalty()] before the
#'   call.
#'
#' @seealso [quadratic_penalty()], [parameters7::kron_identity()]
#'
#' @keywords internal
.kron_quadratic <- function(P, m, link_lambda, tol) {
  e <- eigen(P, symmetric = TRUE)
  keep <- e$values > tol * max(e$values, 0)
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
    logpdet_DPD = m * sum(log(e$values[keep])),
    DPD = Pb,
    tol = tol
  )
}


#' The Quadratic Form of a Quadratic Penalty
#'
#' @description
#' Computes \eqn{(D\beta)'P(D\beta)}, the part of the value that depends on the
#' coefficients. It is shared by [penalty_value.QuadraticPenalty()] and the
#' hyperparameter derivatives, which use \eqn{\beta} only through this
#' quantity.
#'
#' @details
#' Uses \eqn{P} and the map, not the cached \eqn{D'PD}, so the intermediate is
#' the \eqn{m}-vector \eqn{D\beta} and the function does not form a
#' \eqn{q \times q} product.
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
#' Every term is closed form in \eqn{\lambda}, because the rank and the log
#' pseudo-determinant of \eqn{D'PD} are fixed at construction.
#'
#' @details
#' \deqn{\rho(\beta; \lambda) = \frac{\lambda}{2}\,(D\beta)'P(D\beta)
#'   - \frac{r}{2}\log\lambda + \frac{r}{2}\log 2\pi
#'   - \frac{1}{2}\log\mathrm{pdet}(D'PD),}
#'
#' with \eqn{r} the stored rank of \eqn{D'PD}. At \eqn{\beta = 0} the value is
#' the constant alone and still depends on \eqn{\lambda}, so the smoothing
#' parameter remains estimable from a joint objective.
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
#' # At the origin the value is the constant alone, shown at three values of
#' # lambda.
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
    r / 2 * log(2 * pi) - pen@logpdet_DPD / 2
}

#' @title Coefficient Derivatives of a Quadratic Penalty
#' @name penalty_gradient.QuadraticPenalty
#'
#' @description
#' `penalty_gradient()` returns \eqn{\lambda D'PD\beta} and `penalty_hessian()`
#' returns \eqn{\lambda D'PD}. Both read the cached \eqn{D'PD}, so they do not
#' apply the map at call time, and the Hessian does not depend on \eqn{\beta}
#' because the value is a quadratic form.
#'
#' @details
#' Because the Hessian is constant, the gradient is exactly the Hessian applied
#' to the coefficients. When \eqn{D'PD} has full rank, a Newton step on this
#' penalty alone reaches its minimum, at the origin, in one iteration from any
#' starting point. When \eqn{D'PD} is rank deficient, the Hessian is singular
#' and the minimum is attained on the whole null space of \eqn{D'PD}.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A list holding `lambda`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` returns a numeric vector of length `pen@n_coef`;
#'   `penalty_hessian()` a symmetric matrix of that order, a base matrix for an
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
#' The three blocks in \eqn{\lambda} are all closed form. `penalty_grad_theta()`
#' returns \eqn{\tfrac{1}{2}(D\beta)'P(D\beta) - r/(2\lambda)},
#' `penalty_hess_theta()` returns \eqn{r/(2\lambda^2)}, and `penalty_cross()`
#' returns \eqn{D'PD\beta}, which is the gradient with \eqn{\lambda} divided
#' out.
#'
#' @details
#' The value is affine in \eqn{\lambda} apart from the term
#' \eqn{-\tfrac{r}{2}\log\lambda} of the normalizing constant, so the second
#' derivative comes from that term alone and is positive. The value is
#' therefore convex in the smoothing parameter, with its minimum at
#' \eqn{\lambda = r / (D\beta)'P(D\beta)}. That is the joint-mode estimate of
#' \eqn{\lambda} at fixed coefficients, and it is finite because the constant
#' is kept.
#'
#' The mixed block does not contain \eqn{\lambda}, because the coefficient
#' gradient \eqn{\lambda D'PD\beta} is linear in \eqn{\lambda} and its slope
#' \eqn{D'PD\beta} does not depend on \eqn{\lambda}.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A list holding `lambda`.
#' @param scale Read by the generic, which applies the chain rule onto the
#'   unconstrained scale after dispatch. These methods always return the
#'   parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` returns a list of one number named `lambda`;
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
#' This page documents the three predicates that a consumer evaluates before
#' choosing how to fit a quadratic block. `penalty_kinks()` returns
#' `numeric(0)`, because the value is a polynomial and smooth everywhere.
#' `is_quadratic()` returns `TRUE`, so the
#' four marginal quantities are available. `is_proper()` returns `TRUE` only
#' when the penalty has full rank.
#'
#' @details
#' Properness is a rank test and nothing else: \eqn{\exp(-\rho)} is constant
#' along every null direction of \eqn{D'PD}, so it is integrable only when
#' there is no null direction. A ridge over \eqn{q} coefficients is proper;
#' second differences over \eqn{q} coefficients leave the level and the slope
#' unpenalized, have rank \eqn{q - 2}, and give an improper prior. Both can be
#' used, and the difference appears in the normalizing constant, which is
#' taken only over the range.
#'
#' Because there are no kinks, a fit can treat a quadratic block with any
#' method for smooth objectives, and [penalty_prox()] is a linear solve and
#' not a threshold.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param theta A list holding `lambda`. It is accepted by `penalty_kinks()`,
#'   whose result does not depend on it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` returns a numeric vector of length zero.
#'   `is_quadratic()` returns the single logical `TRUE`.
#'   `is_proper()` returns a single logical, `TRUE` when the rank equals the
#'   number of coefficients and the null basis has no columns.
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
#'   generics, [penalty_matrix.QuadraticPenalty()] for the quantities that
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
#' This page documents the four quantities that a REML or marginal likelihood
#' criterion reads. `penalty_matrix()` returns \eqn{\lambda D'PD},
#' `penalty_rank()` returns the rank of \eqn{D'PD} fixed at construction,
#' `penalty_null_basis()` returns the stored orthonormal basis of its null
#' space, and `penalty_logpdet()` returns the log pseudo-determinant with its
#' first two derivatives in \eqn{\lambda}.
#'
#' @details
#' Since \eqn{S(\lambda) = \lambda D'PD} and only the range contributes,
#'
#' \deqn{\log^{+}\lvert S(\lambda) \rvert = r\log\lambda
#'   + \log\mathrm{pdet}(D'PD), \qquad
#'   \frac{\partial}{\partial\lambda} = \frac{r}{\lambda}, \qquad
#'   \frac{\partial^2}{\partial\lambda^2} = -\frac{r}{\lambda^2}.}
#'
#' All three are exact and computed without a decomposition at call time. The
#' rank and the null basis do not depend on \eqn{\lambda} and their generics
#' take no `theta`.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param theta A list holding `lambda`. Only `penalty_matrix()` and
#'   `penalty_logpdet()` read it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_matrix()` returns a symmetric matrix of side `pen@n_coef`,
#'   a base matrix ordinarily and a `dgCMatrix` when the penalty was built with
#'   `blocks > 1`. `penalty_rank()` returns a single integer.
#'   `penalty_null_basis()` returns a matrix with `pen@n_coef` rows and
#'   `pen@n_coef - penalty_rank(pen)` orthonormal columns, sparse under
#'   `blocks > 1`. `penalty_logpdet()` returns a list of `value` (a number),
#'   `grad` (a list of one number, `lambda`) and `hess` (a list of one number,
#'   `lambda_lambda`).
#'
#' @examples
#' pen <- quadratic_penalty(crossprod(diff(diag(5), differences = 2)))
#' penalty_rank(pen)
#'
#' # The null basis spans the straight lines, and the matrix annihilates it.
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
  list(value = r * log(lam) + pen@logpdet_DPD,
       grad = list(lambda = r / lam),
       hess = list(lambda_lambda = -r / lam^2))
}
