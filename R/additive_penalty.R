#' @include quadratic_penalty.R
NULL

#' @title S7 Class for a Sum of Quadratic Penalties
#'
#' @description
#' This class holds the objects that [additive_penalty()] builds: several
#' quadratic penalties added together, each with a smoothing parameter of its
#' own. Beyond the seven properties that every [penalty()] carries, an object
#' of this class stores the list of component matrices, already carried
#' through the map, the rank and null basis of the sum, and the eigenvalue
#' tolerance at which ranks are counted.
#'
#' @details
#' The rank and the null basis are stored because both depend only on the
#' components. The null space of a sum of positive semidefinite
#' matrices with positive weights is the intersection of the null spaces of
#' the components, so the rank and the null basis do not change with the
#' smoothing parameters. A rank counted on the assembled \eqn{S(\lambda)}
#' would fall as the parameters spread apart, because the smaller
#' contributions drop below any fixed relative eigenvalue tolerance; the last
#' example of [additive_penalty()] shows this.
#'
#' @section Methods:
#' The list below gives the methods registered on this class, by the page that
#' documents each.
#'
#' * The value and its coefficient derivatives:
#'   [penalty_value()][penalty_value.AdditivePenalty],
#'   [penalty_gradient()][penalty_gradient.AdditivePenalty] and
#'   [penalty_hessian()][penalty_hessian.AdditivePenalty].
#' * The hyperparameter derivatives:
#'   [penalty_grad_theta()][penalty_grad_theta.AdditivePenalty],
#'   [penalty_hess_theta()][penalty_hess_theta.AdditivePenalty] and
#'   [penalty_cross()][penalty_cross.AdditivePenalty].
#' * The predicates: [penalty_kinks()][penalty_kinks.AdditivePenalty],
#'   [is_proper()][is_proper.AdditivePenalty] and
#'   [is_quadratic()][is_quadratic.AdditivePenalty].
#' * The quantities of a marginal criterion:
#'   [penalty_matrix()][penalty_matrix.AdditivePenalty],
#'   [penalty_rank()][penalty_rank.AdditivePenalty],
#'   [penalty_null_basis()][penalty_null_basis.AdditivePenalty] and
#'   [penalty_logpdet()][penalty_logpdet.AdditivePenalty].
#' * The higher derivatives of a marginal criterion:
#'   [penalty_dhessian()][penalty_dhessian.AdditivePenalty],
#'   [penalty_d2hessian()][penalty_d2hessian.AdditivePenalty],
#'   [penalty_dcross()][penalty_dcross.AdditivePenalty] and
#'   [beta_quadratic()][beta_quadratic.AdditivePenalty].
#' * The movement of the Hessian with the coefficients:
#'   [penalty_dhessian_beta()][penalty_dhessian_beta.AdditivePenalty],
#'   [penalty_d2hessian_beta()][penalty_d2hessian_beta.penalty] and
#'   [penalty_dhessian_beta_theta()][penalty_dhessian_beta_theta.penalty].
#'   The last two return zero matrices, because the Hessian does not depend on
#'   the coefficients, and are documented on the page of the base-class
#'   methods.
#'
#' The class registers no [penalty_prox()] or [penalty_prox_spec()] method,
#' so calls to those generics dispatch to the methods of the base class
#' [penalty()]. It has its own [penalty_draw()] and [penalty_readable()]
#' methods.
#'
#' @inheritParams penalty
#' @param mats A list of symmetric positive semidefinite matrices of the same
#'   side, already carried through the map, so each is `n_coef` by `n_coef`.
#'   One smoothing parameter goes with each, named `lambda1`, `lambda2` and so
#'   on in this order.
#' @param p_rank The rank of the sum, a single whole number, the same at every
#'   positive setting of the parameters.
#' @param null_basis An orthonormal basis of the null space shared by the
#'   components, `n_coef` by `n_coef - p_rank`, read from the same
#'   eigendecomposition as `p_rank` and fixed in the same way.
#' @param tol The relative eigenvalue tolerance, a single positive number,
#'   `1e-10` by default. [additive_penalty()] counts `p_rank` at this
#'   tolerance, and [additive_sum()] passes it to [additive_stable()] as
#'   `r_tol`, so the rank of each group of components is counted at the same
#'   tolerance.
#' @param blocks How many times each component is repeated blockwise, a whole
#'   number. An empty value is read as 1. Above 1, `mats` holds the repeated
#'   components as sparse matrices and `block_mats` holds one block of each.
#' @param block_mats The components of one block, a list of base matrices, or
#'   an empty list when `blocks` is 1.
#' @param readable `"lambda"` or `"sd"`: whether [penalty_readable()] reports
#'   the smoothing parameters themselves or the standard deviations
#'   \eqn{\lambda_k^{-1/2}}. An empty value is read as `"lambda"`.
#'
#' @return An S7 object of class `AdditivePenalty`, inheriting from [penalty()],
#'   with the seven inherited properties and the seven above. Its `map` is
#'   always `NULL`, because a map given to the constructor is absorbed into
#'   `mats`.
#'
#' @seealso [additive_penalty()] for the constructor to use,
#'   [penalty_value.AdditivePenalty()] for the value of the penalty,
#'   [quadratic_penalty()] for one matrix with one scale.
#'
#' @examples
#' pen <- additive_penalty(list(diag(3), diag(c(1, 0, 0))))
#' S7::S7_inherits(pen, AdditivePenalty)
#'
#' # The penalty has one hyperparameter per component and a rank fixed at
#' # construction.
#' pen@params
#' c(rank = pen@p_rank, null_width = ncol(pen@null_basis))
#' length(pen@mats)
#'
#' @keywords internal
#' @export
AdditivePenalty <- S7::new_class(
  name = "AdditivePenalty",
  parent = penalty,
  properties = list(
    mats = S7::class_list,
    p_rank = S7::class_numeric,
    null_basis = S7::class_any,
    tol = S7::new_property(S7::class_numeric, default = 1e-10),
    blocks = S7::class_numeric,
    block_mats = S7::class_list,
    readable = S7::class_character
  )
)

#' How Many Blocks an Additive Penalty Repeats
#'
#' @description
#' Returns the `blocks` property as a single integer, reading an empty value
#' as 1.
#'
#' @param pen An [AdditivePenalty()] object.
#'
#' @return A single integer of at least 1.
#'
#' @keywords internal
additive_blocks <- function(pen) {
  b <- pen@blocks
  if (!length(b)) 1L else as.integer(b[[1L]])
}

#' @title Construct a Sum of Quadratic Penalties
#'
#' @description
#' Builds the Gaussian prior whose precision is a weighted sum of fixed
#' positive semidefinite matrices, with one smoothing parameter for each
#' matrix. An anisotropic tensor-product smooth uses this form: with one
#' parameter per margin, the fit can be rough in one direction and smooth in
#' another, whereas a single parameter smooths every direction alike.
#'
#' @details
#' # The value
#'
#' \deqn{\rho(\beta;\lambda) = \tfrac{1}{2}\beta^\top S(\lambda)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,S(\lambda)
#'   + \tfrac{r}{2}\log 2\pi, \qquad
#'   S(\lambda) = \sum_{k} \lambda_k P_k,}
#'
#' with \eqn{r} the rank of the sum. The hyperparameters are named `lambda1`,
#' `lambda2` and so on, in the order in which the components were given. Each
#' is positive and is carried onto the real line by `link_lambda`, a log link
#' by default.
#'
#' # The derivatives of the log pseudo-determinant
#'
#' With \eqn{S^{+}} the pseudo-inverse of \eqn{S(\lambda)}, the derivatives
#' that a marginal criterion reads are traces against it:
#'
#' \deqn{\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\,S
#'   = \operatorname{tr}(S^{+}P_k), \qquad
#'   \frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
#'   \log\mathrm{pdet}\,S = -\operatorname{tr}(S^{+}P_kS^{+}P_l).}
#'
#' Both formulas are exact. [additive_sum()] evaluates the traces, together
#' with the log pseudo-determinant, by the similarity transformation of
#' [additive_stable()], and never forms \eqn{S^{+}} itself.
#'
#' # The rank of the sum
#'
#' The rank and the null basis are computed once, at construction, from the
#' sum of the components after each has been divided by its largest absolute
#' entry, and the rank is the number of eigenvalues of that sum above `tol`
#' times the largest. Both are then fixed: the null space of \eqn{S(\lambda)} is
#' the intersection of the null spaces of the components at every positive
#' \eqn{\lambda}. An eigenvalue count of the assembled \eqn{S(\lambda)} itself
#' is correct only while the parameters are comparable. Once they differ by many
#' orders of magnitude, the smaller contributions fall below any fixed relative
#' tolerance and the count drops, although the stored null basis still
#' annihilates the sum to rounding error.
#'
#' # Parameters that differ by many orders of magnitude
#'
#' The stored rank fixes how many eigenvalues enter the log
#' pseudo-determinant, but it does not make those eigenvalues computable from
#' the assembled sum. When the parameters differ by enough orders of magnitude,
#' the condition number of \eqn{S(\lambda)} exceeds \eqn{1/\epsilon} and the
#' smallest eigenvalue on the range falls below the absolute accuracy of one
#' symmetric eigendecomposition, which is of order
#' \eqn{d\,\epsilon\lVert S\rVert_{2}}. The log pseudo-determinant is therefore
#' not taken, in general, from one decomposition of the assembled sum.
#' [additive_stable()], which [additive_sum()] calls, partitions the
#' components by the size of their contributions and transforms by
#' similarity, so that each decomposition resolves only the spread within one
#' group of components. The sum is decomposed as a whole only when every
#' contribution is at least [additive_tol()] times the largest, or when the
#' dominant group by itself has the rank of the sum. The identities of the
#' transformation are exact, and the accuracy of the value
#' and of both derivative blocks is limited by the conditioning of the
#' dominant group at each level of the partition rather than by the spread of
#' the parameters. Where the weighted sum or one of its reductions overflows,
#' or where the rank identity that [additive_sum()] checks at every call is
#' violated or cannot be evaluated, the log pseudo-determinant and its
#' derivatives are `NaN`, without an error. Such points become more frequent
#' as the spread and the number of components grow.
#'
#' # The quantities of a marginal criterion
#'
#' [is_quadratic()] returns `TRUE`, because a sum of quadratic forms is a
#' quadratic form, so [penalty_matrix()], [penalty_rank()],
#' [penalty_null_basis()] and [penalty_logpdet()] all have methods for this
#' class, and [check_penalty()] runs its quadratic rows here (the row of the
#' null basis only when the null basis has columns). The matrix moves with one
#' hyperparameter per component, whereas the plain quadratic branch has a single
#' scale, so the log pseudo-determinant is not linear in the logarithm of any
#' one parameter. The row of [check_penalty()] that checks it therefore compares
#' its gradient with `numDeriv`, where the plain quadratic branch checks a
#' slope.
#'
#' [has_prox()] returns `FALSE`. That function first tests whether the class
#' registers a [penalty_prox()] method of its own, before it calls
#' [is_quadratic()], and this class registers none.
#'
# Blockwise repetition
#'
#' `blocks = m` builds the penalty whose components are
#' \eqn{I_m \otimes P_k}, which is the form of a smooth repeated once per
#' group. The sum is \eqn{I_m \otimes S_1(\lambda)}, with \eqn{S_1} the sum
#' over one block, so its rank is \eqn{m} times the rank of the block, its log
#' pseudo-determinant is \eqn{m\log\mathrm{pdet}\,S_1(\lambda)}, and both
#' derivative blocks are \eqn{m} times those of the block. The decomposition is
#' taken on one block, and the repeated components are stored as sparse
#' matrices, which is also the storage that [penalty_hessian()] and
#' [penalty_matrix()] return. `blocks` cannot be combined with `map`.
#'
#' # The readable quantities
#'
#' Setting `readable = "sd"` declares the smoothing parameters to be
#' precisions of a prior, and [penalty_readable()] then reports the standard
#' deviations \eqn{\sigma_k = \lambda_k^{-1/2}}, named `sigma1`, `sigma2` and
#' so on, in place of the parameters. The penalty itself is unchanged.
#'
#' @param mats A list of at least one symmetric positive semidefinite matrix,
#'   all with the same number of rows. Each component is symmetrized. A
#'   component is rejected with an error that gives its position when it has
#'   an entry that is not finite (`NA`, `NaN` or infinite), when it differs
#'   from its transpose by more than `1e-8` times the larger of 1 and its
#'   largest absolute entry, when its largest eigenvalue is not positive (a
#'   zero or negative semidefinite matrix), or when its smallest eigenvalue is
#'   below `-tol` times its largest. Under `blocks > 1` they are the
#'   components of one block.
#' @param map The matrix \eqn{D}, with as many rows as each component has and
#'   one column per coefficient, or `NULL` (the default) for the identity. A
#'   map is absorbed at construction, and each component becomes
#'   \eqn{D^\top P_kD}, so the stored `map` is always `NULL` and `n_coef` is
#'   the number of columns of the map. A map with an entry that is not finite
#'   is rejected, and so is a map that carries a component to the zero matrix,
#'   with an error that gives the position of that component.
#' @param link_lambda The \pkg{linkfunctions7} link that carries each smoothing
#'   parameter onto the whole real line, `linkfunctions7::log_link()` by
#'   default. The same link is used for every parameter.
#' @param blocks How many times each component is repeated blockwise, a whole
#'   number of at least 1. `1L`, the default, gives the ordinary penalty. A
#'   value above 1 cannot be combined with `map`.
#' @param readable `"lambda"` (the default) or `"sd"`, the quantities that
#'   [penalty_readable()] reports. See the section on the readable
#'   quantities.
#' @param tol The relative eigenvalue tolerance, `1e-10` by default. It is used
#'   to reject a component whose smallest eigenvalue is below `-tol` times its
#'   largest, and to count the rank of the sum of the components after each has
#'   been divided by its largest absolute entry. The object stores `tol`, and
#'   [additive_sum()] passes it to [additive_stable()] as `r_tol`, so the rank
#'   of each group of components in the similarity transformation is counted
#'   at the same tolerance as the rank of the sum.
#'
#' @return An [AdditivePenalty()] object with one hyperparameter per component,
#'   named `lambda1`, `lambda2`, ... The hyperparameters are bounded on
#'   \eqn{(0, \infty)}.
#'
#' @examples
#' # This is an anisotropic tensor smooth over a 4 x 4 grid, and the curvature
#' # along each margin is penalized separately.
#' P <- crossprod(diff(diag(4), differences = 2))
#' pen <- additive_penalty(list(kronecker(diag(4), P), kronecker(P, diag(4))))
#' pen@params
#' c(n_coef = pen@n_coef, rank = penalty_rank(pen))
#'
#' # Raising one parameter penalizes roughness along that margin alone.
#' set.seed(2)
#' b <- rnorm(16)
#' penalty_value(pen, b, list(lambda1 = 1, lambda2 = 1))
#' penalty_value(pen, b, list(lambda1 = 1, lambda2 = 100))
#'
#' # With one component the value is that of the plain quadratic penalty.
#' P2 <- crossprod(diff(diag(5), differences = 2))
#' bb <- c(0.4, -1.1, 0.7, 0.2, -0.3)
#' penalty_value(additive_penalty(list(P2)), bb, list(lambda1 = 3)) -
#'   penalty_value(quadratic_penalty(P2), bb, list(lambda = 3))
#'
#' # The stored rank stays fixed as the parameters spread apart, while an
#' # eigenvalue count of the assembled sum drops.
#' a <- additive_penalty(list(crossprod(diff(diag(5))), P2))
#' penalty_rank(a)
#' sapply(c(1, 1e10), function(r) {
#'   S <- penalty_matrix(a, list(lambda1 = 1, lambda2 = r))
#'   ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
#'   sum(ev > 1e-10 * max(ev))
#' })
#'
#' # blocks = m repeats each component without decomposing the repetition.
#' rep4 <- additive_penalty(list(crossprod(diff(diag(5))), P2), blocks = 4)
#' c(n_coef = rep4@n_coef, rank = penalty_rank(rep4))
#' th <- list(lambda1 = 2, lambda2 = 0.5)
#' penalty_logpdet(rep4, th)$value - 4 * penalty_logpdet(a, th)$value
#'
#' # readable = "sd" reports the standard deviations of the prior.
#' sdp <- additive_penalty(list(diag(2), diag(c(1, 2))), readable = "sd")
#' penalty_readable(sdp, list(lambda1 = 4, lambda2 = 0.25))$value
#'
#' @seealso [quadratic_penalty()] for one matrix with one scale,
#'   [structured_penalty()] for a matrix whose entries move with the
#'   hyperparameters, [penalty_logpdet.AdditivePenalty()] for the log
#'   pseudo-determinant and its derivatives, [additive_sum()] for the
#'   computation that they share.
#'
#' @export
additive_penalty <- function(mats, map = NULL, blocks = 1L,
                             readable = c("lambda", "sd"),
                             link_lambda = linkfunctions7::log_link(),
                             tol = 1e-10) {
  if (!is.list(mats) || length(mats) < 1L) {
    stop("'mats' must be a non-empty list of matrices.", call. = FALSE)
  }
  readable <- match.arg(readable)
  if (!is.numeric(blocks) || length(blocks) != 1L || is.na(blocks) ||
      blocks < 1 || blocks != round(blocks)) {
    stop("'blocks' must be a whole number of at least 1.", call. = FALSE)
  }
  blocks <- as.integer(blocks)
  if (blocks > 1L && !is.null(map)) {
    stop(paste0("'blocks' and 'map' do not combine: a map mixes the blocks,",
                " so the\n  components are not block diagonal with the",
                " structure that 'blocks' names."), call. = FALSE)
  }
  mats <- lapply(seq_along(mats), function(k) {
    P <- as.matrix(mats[[k]])
    if (!all(is.finite(P))) {
      stop(sprintf("component %d of 'mats' must have only finite entries.", k),
           call. = FALSE)
    }
    if (nrow(P) != ncol(P) ||
        max(abs(P - t(P))) > 1e-8 * max(1, max(abs(P)))) {
      stop(sprintf("component %d of 'mats' must be a symmetric matrix.", k),
           call. = FALSE)
    }
    (P + t(P)) / 2
  })
  d <- unique(vapply(mats, nrow, integer(1)))
  if (length(d) != 1L) {
    stop("every component of 'mats' must have the same dimension.",
         call. = FALSE)
  }
  for (k in seq_along(mats)) {
    ev <- eigen(mats[[k]], symmetric = TRUE, only.values = TRUE)$values
    if (max(ev) <= 0) {
      stop(sprintf("component %d is zero or negative semidefinite.", k),
           call. = FALSE)
    }
    if (min(ev) < -tol * max(ev)) {
      stop(sprintf("component %d is not positive semidefinite.", k),
           call. = FALSE)
    }
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
    if (nrow(map) != d) {
      stop(sprintf("'map' must have %d rows.", d), call. = FALSE)
    }
    # densified for the reason quadratic_penalty() states: a Matrix map would
    # otherwise decide the class of every component and of every quantity
    # assembled from them
    mats <- lapply(mats, function(P) as.matrix(crossprod(map, P %*% map)))
    d <- ncol(map)
    # a map can annihilate a component that was nonzero before it, and a zero
    # component cannot be normalized for the rank below
    for (k in seq_along(mats)) {
      ev <- eigen(mats[[k]], symmetric = TRUE, only.values = TRUE)$values
      if (max(ev) <= 0) {
        stop(sprintf("component %d is zero after 'map' is applied.", k),
             call. = FALSE)
      }
    }
  }

  # The rank is a property of the components, not of any one value of the
  # parameters: the null space of the sum is the intersection of the null
  # spaces, so the normalized components are added once and its rank read
  # there. Reading it off S(lambda) would make the rank fall as the
  # parameters spread apart.
  stacked <- Reduce(`+`, lapply(mats, function(P) P / max(abs(P))))
  es <- eigen(stacked, symmetric = TRUE)
  keep <- es$values > tol * max(es$values)
  r <- sum(keep)
  # the null space of the sum is the intersection of the components', so
  # these vectors span the null space of S(lambda) at every admissible
  # lambda, which is what makes it a property of the family
  nb <- es$vectors[, !keep, drop = FALSE]

  # REPEATED BLOCKWISE, the rank, the null basis and every determinant follow
  # from one block, the eigenvalues of I_m (x) S_1 being S_1's repeated m
  # times; the repeated components are stored sparse and never decomposed
  bm <- list()
  if (blocks > 1L) {
    bm <- mats
    I <- Matrix::Diagonal(blocks)
    mats <- lapply(bm, function(P) Matrix::kronecker(
      I, methods::as(Matrix::Matrix(P, sparse = TRUE), "generalMatrix")))
    nb <- Matrix::kronecker(I, Matrix::Matrix(nb, sparse = TRUE))
    r <- blocks * r
    d <- blocks * d
  }

  nm <- paste0("lambda", seq_along(mats))
  AdditivePenalty(
    penalty_name = sprintf("additive [%d components]", length(mats)),
    map = NULL,
    n_coef = d,
    params = nm,
    params_bounds = stats::setNames(rep(list(c(0, Inf)), length(nm)), nm),
    link_params = stats::setNames(
      replicate(length(nm), link_lambda, simplify = FALSE), nm),
    params_smooth = stats::setNames(rep(TRUE, length(nm)), nm),
    mats = mats,
    p_rank = r,
    null_basis = nb,
    tol = tol,
    blocks = blocks,
    block_mats = bm,
    readable = readable
  )
}

#' The Weighted Sum and Its Log Pseudo-Determinant
#'
#' @description
#' Assembles \eqn{S(\lambda) = \sum_k \lambda_k P_k} and returns it together
#' with its log pseudo-determinant and the two derivative blocks
#' \eqn{\operatorname{tr}(S^{+}P_k)} and
#' \eqn{\operatorname{tr}(S^{+}P_kS^{+}P_l)}. The value, the hyperparameter
#' derivatives and [penalty_logpdet()] of an [AdditivePenalty()] read these
#' quantities from this function, which calls [additive_stable()] for the
#' last three.
#'
#' @details
#' # A single decomposition of the sum
#'
#' A single eigendecomposition of \eqn{S(\lambda)}, keeping its `p_rank`
#' largest eigenvalues, fails once the parameters spread apart. A symmetric
#' eigendecomposition computes its eigenvalues with an absolute accuracy of
#' order \eqn{d\,\epsilon\lVert S\rVert_{2}}, so the smallest eigenvalue on the
#' range can fall below that resolution and carry no significant digit. The
#' sum of logarithms is then `NaN` if that eigenvalue is computed as negative
#' and an inaccurate number if it is computed as positive, and the computation
#' raises no warning in either case.
#'
#' # The similarity transformation
#'
#' [additive_stable()] partitions the components by the size of their
#' contributions, \eqn{\lambda_k\max_{ij}\lvert(P_k)_{ij}\rvert}, and rotates
#' onto the eigenvectors of the dominant group, the components whose
#' contribution is at least [additive_tol()] times the largest. Let
#' \eqn{U_{+}} span the range of that group and \eqn{U_{0}} its kernel, and let
#' \eqn{A = U_{+}^{\top}SU_{+}}, \eqn{C = U_{+}^{\top}SU_{0}},
#' \eqn{F = A^{-1}C} and \eqn{M = U_{0}^{\top}SU_{0} - C^{\top}F}. Then
#'
#' \deqn{\log\mathrm{pdet}\,S = \log\lvert A\rvert + \log\mathrm{pdet}\,M,
#'   \qquad S^{+} = U_{+}A^{-1}U_{+}^{\top} + GM^{+}G^{\top},
#'   \qquad G = U_{0} - U_{+}F.}
#'
#' Both identities are exact, because the null space of \eqn{S} lies in the
#' span of \eqn{U_{0}}. In the first, \eqn{A} is the whole sum restricted to
#' the range of the dominant group, where the dominant group sets the scale,
#' and \eqn{M} is a Schur complement of the size of the subordinate
#' components, so the decomposition of \eqn{A} does not have to resolve the
#' spread between the two. The second writes every trace against \eqn{S^{+}}
#' as a trace against \eqn{A^{-1}} plus a trace of the same kind against
#' \eqn{M^{+}}, so the derivatives are taken in the transformed coordinates
#' and never from a materialized \eqn{S^{+}}. \eqn{M} is treated in the same
#' way at the next level. The recursion stops at a level where every
#' remaining component joins the dominant group, or where the rank of the
#' dominant group is zero or reaches the rank left, and that level takes one
#' eigendecomposition.
#'
#' Two quantities are exactly zero and are dropped instead of computed. The
#' dominant group vanishes on \eqn{U_{0}}, because that subspace is its kernel,
#' so a computed value there would consist of rounding errors of order
#' \eqn{\lambda_{\mathrm{dom}}\epsilon}, which would swamp the subordinate
#' terms; the blocks of each dominant component on \eqn{U_{0}} vanish for the
#' same reason. The part of \eqn{M} that the next level does not carry as
#' components is likewise built from the reductions of the dominant
#' components and not obtained by subtracting the subordinate reductions
#' from \eqn{M}, since that subtraction would be a difference of two
#' quantities of the dominant size.
#'
#' # The rank identity
#'
#' \eqn{\operatorname{tr}(S^{+}S) = r} exactly, so
#' \eqn{\sum_k\lambda_k\operatorname{tr}(S^{+}P_k)} must equal the rank. The
#' identity needs no reference and is checked at every call. Where the value
#' returned by [additive_stable()] is not finite, where the identity is
#' violated by more than `additive_check_tol()` times the rank, or where it
#' cannot be evaluated because a trace overflows, the point is treated as not
#' resolvable in double precision, and `logpdet`, `dlog` and `d2log` are
#' `NaN`, without an error. `S` is assembled without any decomposition and is
#' returned as it is. [penalty_matrix()], [penalty_gradient()] and
#' [penalty_hessian()] assemble the sum themselves and do not call this
#' function, so they stay finite at such a point unless the weighted sum
#' itself overflows.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param theta The aligned hyperparameter list, as [align_ptheta()] returns
#'   it.
#'
#' @return A list of four: `S`, the assembled symmetric matrix of side
#'   `pen@n_coef`; `logpdet`, a single number; `dlog`, a numeric vector of
#'   \eqn{\operatorname{tr}(S^{+}P_k)} in `pen@params` order; and `d2log`, the
#'   square matrix of \eqn{\operatorname{tr}(S^{+}P_kS^{+}P_l)}. The last
#'   three are `NaN` at a setting that the transformation cannot resolve.
#'
#' @seealso [additive_stable()] for the transformation,
#'   [additive_penalty()], [penalty_logpdet.AdditivePenalty()]
#'
#' @keywords internal
additive_sum <- function(pen, theta) {
  lam <- unlist(theta[pen@params])
  S <- Reduce(`+`, Map(function(P, l) l * P, pen@mats, lam))
  m <- additive_blocks(pen)
  st <- if (m > 1L) {
    # one block decomposed, and every determinant quantity m times its own
    s1 <- additive_stable(pen@block_mats, lam, pen@p_rank / m,
                          r_tol = pen@tol)
    list(value = m * s1$value, dlog = m * s1$dlog, d2log = m * s1$d2log)
  } else {
    additive_stable(pen@mats, lam, pen@p_rank, r_tol = pen@tol)
  }
  # tr(S+ S) = r exactly, so this sum is the rank whatever the parameters are.
  # It is free -- the gradient is computed anyway -- it needs no reference,
  # and it is the only check available where the spread puts a finite
  # difference of the value out of reach.
  # An overflowing trace makes the comparison NA, which is as unresolved as
  # a violated identity.
  good <- is.finite(st$value) &&
    abs(sum(lam * st$dlog) - pen@p_rank) <=
      additive_check_tol() * pen@p_rank
  if (!isTRUE(good)) {
    K <- length(lam)
    return(list(S = S, logpdet = NaN, dlog = rep(NaN, K),
                d2log = matrix(NaN, K, K)))
  }
  list(S = S, logpdet = st$value, dlog = st$dlog, d2log = st$d2log)
}

#' The Tolerance of the Rank Identity
#'
#' @description
#' Returns `1e-3`, the relative departure from
#' \eqn{\sum_k\lambda_k\operatorname{tr}(S^{+}P_k) = r} beyond which
#' [additive_sum()] reports `NaN` instead of a number. The departure is
#' taken relative to the rank \eqn{r}.
#'
#' @details
#' The identity is exact, so a departure from it is accumulated rounding in
#' the similarity transformation of [additive_stable()]. The threshold is
#' chosen above the departures that rounding produces at points that the
#' transformation resolves, so that a departure beyond it signals that the
#' smaller eigenvalues have been lost. The threshold is not a bound on the
#' error of the returned quantities: at a point that passes the check, the
#' error of the log pseudo-determinant can exceed the departure from the
#' identity.
#'
#' @return A single number.
#'
#' @seealso [additive_sum()], which reads it.
#'
#' @keywords internal
additive_check_tol <- function() 1e-3

#' The Split Threshold of the Components
#'
#' @description
#' Returns \eqn{\epsilon^{0.4}} (about `5.5e-07`), with \eqn{\epsilon} the
#' machine precision. In [additive_stable()], a component whose contribution
#' \eqn{\lambda_k\max_{ij}\lvert(P_k)_{ij}\rvert} is at least this fraction of
#' the largest contribution joins the dominant group, and the other components
#' are carried to the next level.
#'
#' @details
#' The threshold balances two sources of rounding. A smaller value widens the
#' dominant group, so one decomposition must resolve a spread of up to the
#' reciprocal of the threshold. A larger value narrows the group, so the
#' recursion runs deeper and every level adds the rounding of one more
#' reduction. The exponent 0.4 was chosen on the rank identity of
#' [additive_check_tol()]; \pkg{mgcv}'s `gam.reparam` uses
#' \eqn{\epsilon^{0.3}} in the same role.
#'
#' @return A single number.
#'
#' @seealso [additive_stable()], which reads it.
#'
#' @keywords internal
additive_tol <- function() .Machine$double.eps^0.4

#' The Stable Log Pseudo-Determinant and Its Derivative Blocks
#'
#' @description
#' Evaluates \eqn{\log\mathrm{pdet}\sum_k\lambda_kP_k} together with
#' \eqn{\operatorname{tr}(S^{+}P_k)} and
#' \eqn{\operatorname{tr}(S^{+}P_kS^{+}P_l)} by the similarity transformation
#' that [additive_sum()] describes. The function recurses on the Schur
#' complement, so it handles any number of components.
#'
#' @details
#' The rank of a group of components is counted on the sum of its
#' components after each has been divided by its largest absolute entry, and
#' never on the weighted sum. [additive_penalty()] counts the rank of the whole
#' penalty in the same way and at the same tolerance. The null space of a sum
#' of positive semidefinite matrices is the intersection of their null spaces,
#' so the rank of a group is a property of that group and is counted without
#' the spread of the parameters. The rank passed to the next level is the
#' current rank less the rank of the dominant group, so it is obtained by
#' subtraction and never counted again.
#'
#' A level returns `NaN` for the value and for both derivative blocks when its
#' weighted sum or reduction is not finite, when the Cholesky factorization of
#' \eqn{A} fails, or, at the last level, when the smallest of the eigenvalues
#' that are kept is not positive; a `NaN` at a deeper level is passed up.
#'
#' @param mats The component matrices, as [AdditivePenalty()] stores them.
#' @param lam The smoothing parameters, in the same order.
#' @param p_rank The rank of the sum, fixed at construction.
#' @param d_tol The threshold at which the components are split,
#'   [additive_tol()] by default.
#' @param r_tol The relative eigenvalue tolerance at which the rank of a group
#'   is counted, `1e-10` by default. [additive_sum()] passes the `tol` stored
#'   in the [AdditivePenalty()] object, which is the tolerance at which
#'   [additive_penalty()] counted the rank of the sum.
#'
#' @return A list of `value` (a single number), `dlog` (a numeric vector),
#'   `d2log` (a square matrix) and `depth` (the number of levels of the
#'   partition, a diagnostic).
#'
#' @seealso [additive_sum()], which calls it and checks its result.
#'
#' @keywords internal
additive_stable <- function(mats, lam, p_rank, d_tol = additive_tol(),
                            r_tol = 1e-10) {
  K <- length(mats)
  grp_rank <- function(Q) {
    St <- Reduce(`+`, lapply(Q, function(P) P / max(abs(P))))
    ev <- eigen(St, symmetric = TRUE)$values
    sum(ev > r_tol * max(ev))
  }
  sym <- function(Z) (Z + t(Z)) / 2

  rec <- function(R, peeled, E, r) {
    act <- which(!peeled)
    cur <- Reduce(`+`, Map(function(P, w) w * P, R[act], lam[act]))
    if (!is.null(E)) cur <- cur + E
    cur <- sym(cur)
    # a sum or a reduction that overflowed cannot be decomposed
    if (!all(is.finite(cur))) {
      return(list(value = NaN, g = rep(NaN, K), T2 = matrix(NaN, K, K),
                  depth = 1L, tr1 = function(X) NaN))
    }
    scal <- lam[act] * vapply(R[act], function(P) max(abs(P)), numeric(1))
    dom <- act[scal >= d_tol * max(scal)]
    ra <- if (length(dom) == length(act)) r else min(grp_rank(R[dom]), r)

    if (length(dom) == length(act) || ra == 0L || ra >= r) {
      e <- eigen(cur, symmetric = TRUE)
      k <- order(e$values, decreasing = TRUE)[seq_len(r)]
      if (e$values[k[r]] <= 0) {
        return(list(value = NaN, g = rep(NaN, K),
                    T2 = matrix(NaN, K, K), depth = 1L,
                    tr1 = function(X) NaN))
      }
      V <- e$vectors[, k, drop = FALSE]
      Sp <- V %*% (t(V) / e$values[k])
      SpR <- lapply(R, function(P) Sp %*% P)
      T2 <- matrix(0, K, K)
      for (i in seq_len(K)) for (j in i:K) {
        T2[i, j] <- T2[j, i] <- sum(t(SpR[[i]]) * SpR[[j]])
      }
      return(list(value = sum(log(e$values[k])), depth = 1L,
                  g = vapply(R, function(P) sum(Sp * P), numeric(1)),
                  T2 = T2, tr1 = function(X) sum(Sp * X)))
    }

    Sa <- Reduce(`+`, Map(function(P, w) w * P, R[dom], lam[dom]))
    ea <- eigen(sym(Sa), symmetric = TRUE)
    o <- order(ea$values, decreasing = TRUE)
    Up <- ea$vectors[, o[seq_len(ra)], drop = FALSE]
    U0 <- ea$vectors[, o[-seq_len(ra)], drop = FALSE]

    # the dominant group is EXACTLY zero on U0, that subspace being its
    # kernel, so the subordinate components alone are read there
    sub <- setdiff(act, dom)
    subm <- Reduce(`+`, Map(function(P, w) w * P, R[sub], lam[sub]))
    if (!is.null(E)) subm <- subm + E
    A <- sym(crossprod(Up, cur %*% Up))
    Cm <- crossprod(Up, subm %*% U0)
    cA <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(cA)) {
      return(list(value = NaN, g = rep(NaN, K), T2 = matrix(NaN, K, K),
                  depth = 1L, tr1 = function(X) NaN))
    }
    Ainv <- chol2inv(cA)
    Fm <- Ainv %*% Cm

    ak <- lapply(R, function(P) crossprod(Up, P %*% Up))
    red <- lapply(seq_len(K), function(k) {
      if (k %in% dom) {
        # its own U0 blocks are zero for the same reason: dropped
        sym(crossprod(Fm, ak[[k]] %*% Fm))
      } else {
        ck <- crossprod(Up, R[[k]] %*% U0)
        bk <- crossprod(U0, R[[k]] %*% U0)
        sym(bk - crossprod(ck, Fm) - crossprod(Fm, ck) +
              crossprod(Fm, ak[[k]] %*% Fm))
      }
    })
    # BUILT, not subtracted: M is the sum of every active component's
    # reduction plus G'EG, so what the next level does not carry is the
    # dominant reductions plus G'EG, each of the subordinate size.
    E2 <- Reduce(`+`, Map(function(P, w) w * P, red[dom], lam[dom]))
    if (!is.null(E)) {
      Ea <- crossprod(Up, E %*% Up)
      Ec <- crossprod(Up, E %*% U0)
      Eb <- crossprod(U0, E %*% U0)
      E2 <- E2 + sym(Eb - crossprod(Ec, Fm) - crossprod(Fm, Ec) +
                       crossprod(Fm, Ea %*% Fm))
    }
    peel2 <- peeled
    peel2[dom] <- TRUE
    nx <- rec(red, peel2, E2, r - ra)
    if (!is.finite(nx$value)) {
      return(list(value = NaN, g = rep(NaN, K), T2 = matrix(NaN, K, K),
                  depth = nx$depth + 1L, tr1 = function(X) NaN))
    }

    JX <- lapply(seq_len(K), function(k) {
      if (k %in% dom) -(ak[[k]] %*% Fm)
      else crossprod(Up, R[[k]] %*% U0) - ak[[k]] %*% Fm
    })
    AiA <- lapply(ak, function(Z) Ainv %*% Z)
    T2 <- matrix(0, K, K)
    for (i in seq_len(K)) for (j in i:K) {
      mid <- crossprod(JX[[j]], Ainv %*% JX[[i]])
      T2[i, j] <- T2[j, i] <- sum(t(AiA[[i]]) * AiA[[j]]) +
        nx$tr1(mid + t(mid)) + nx$T2[i, j]
    }
    list(
      value = 2 * sum(log(diag(cA))) + nx$value, depth = nx$depth + 1L,
      g = vapply(seq_len(K), function(k) sum(Ainv * ak[[k]]) + nx$g[k],
                 numeric(1)),
      T2 = T2,
      tr1 = function(X) {
        Xa <- crossprod(Up, X %*% Up)
        Xc <- crossprod(Up, X %*% U0)
        Xb <- crossprod(U0, X %*% U0)
        sum(Ainv * Xa) + nx$tr1(sym(Xb - crossprod(Xc, Fm) -
                                      crossprod(Fm, Xc) +
                                      crossprod(Fm, Xa %*% Fm)))
      })
  }

  st <- rec(mats, rep(FALSE, K), NULL, p_rank)
  list(value = st$value, dlog = st$g, d2log = st$T2, depth = st$depth)
}

#' @title Value of an Additive Penalty
#' @name penalty_value.AdditivePenalty
#'
#' @description
#' Returns the negative log-density of the Gaussian prior whose precision is
#' the weighted sum \eqn{S(\lambda) = \sum_k \lambda_k P_k}. The normalizing
#' constant is taken over the range of the sum.
#'
#' @details
#' \deqn{\rho(\beta;\lambda) = \tfrac{1}{2}\beta^\top S(\lambda)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,S(\lambda)
#'   + \tfrac{r}{2}\log 2\pi,}
#'
#' with \eqn{r} the stored rank. [additive_sum()] assembles the matrix
#' \eqn{S(\lambda)} of the quadratic form and evaluates the log
#' pseudo-determinant by a similarity transformation, which decomposes the
#' assembled sum as a whole only when every component joins the dominant
#' group or when the dominant group by itself has the rank of the sum. The
#' accuracy of the log pseudo-determinant is therefore limited by
#' the conditioning of the dominant group of components at each level of the
#' transformation, and not by the spread of the parameters.
#'
#' With a single component the value equals that of [quadratic_penalty()] up
#' to rounding error.
#'
#' The value is `NaN`, without an error, at a setting that [additive_sum()]
#' reports as not resolvable: where the weighted sum or one of its reductions
#' overflows, or where the rank identity of [additive_check_tol()] is violated
#' or cannot be evaluated.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A named list of `lambda1`, `lambda2`, ..., already aligned and
#'   bound-checked by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number.
#'
#' @examples
#' P1 <- crossprod(diff(diag(5)))
#' P2 <- crossprod(diff(diag(5), differences = 2))
#' pen <- additive_penalty(list(P1, P2))
#' b <- c(0.4, -1.1, 0.7, 0.2, -0.3)
#'
#' penalty_value(pen, b, list(lambda1 = 2, lambda2 = 0.5))
#'
#' # A constant vector lies in the null space of both components, so only the
#' # normalizing constant remains.
#' penalty_value(pen, rep(1, 5), list(lambda1 = 2, lambda2 = 0.5))
#' penalty_value(pen, rep(0, 5), list(lambda1 = 2, lambda2 = 0.5))
#'
#' @seealso [additive_penalty()] for the construction,
#'   [penalty_gradient.AdditivePenalty()] for the coefficient derivatives,
#'   [additive_sum()] for the quantities that it reads.
#' @keywords internal
S7::method(penalty_value, AdditivePenalty) <- function(pen, beta, theta, ...) {
  a <- additive_sum(pen, theta)
  0.5 * sum(beta * as.numeric(a$S %*% beta)) - 0.5 * a$logpdet +
    pen@p_rank / 2 * log(2 * pi)
}

#' @title Coefficient Derivatives of an Additive Penalty
#' @name penalty_gradient.AdditivePenalty
#'
#' @description
#' `penalty_gradient()` returns \eqn{S(\lambda)\beta} and `penalty_hessian()`
#' returns \eqn{S(\lambda)}. Both assemble the weighted sum directly and take
#' no eigendecomposition, because the normalizing constant does not depend on
#' the coefficients.
#'
#' @details
#' The value is a quadratic form, so the Hessian does not move with the
#' coefficients and the gradient is the Hessian applied to them. A coefficient
#' vector in the intersection of the null spaces of the components has zero
#' gradient at every setting of the parameters.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of `lambda1`, `lambda2`, ...
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` returns a numeric vector of length `pen@n_coef`;
#'   `penalty_hessian()` a symmetric base matrix of that side.
#'
#' @examples
#' P1 <- crossprod(diff(diag(5)))
#' P2 <- crossprod(diff(diag(5), differences = 2))
#' pen <- additive_penalty(list(P1, P2))
#' th <- list(lambda1 = 2, lambda2 = 0.5)
#' b <- c(0.4, -1.1, 0.7, 0.2, -0.3)
#'
#' all.equal(penalty_gradient(pen, b, th),
#'           drop(penalty_hessian(pen, b, th) %*% b))
#'
#' # The constants are in both null spaces, so they cost nothing at any
#' # setting of the two parameters.
#' penalty_gradient(pen, rep(1, 5), th)
#' penalty_gradient(pen, rep(1, 5), list(lambda1 = 1e6, lambda2 = 1e-6))
#'
#' @seealso [penalty_value.AdditivePenalty()] for the quantity differentiated,
#'   [penalty_grad_theta.AdditivePenalty()] for the hyperparameter blocks.
#' @keywords internal
S7::method(penalty_gradient, AdditivePenalty) <- function(pen, beta, theta, ...) {
  lam <- unlist(theta[pen@params])
  as.numeric(Reduce(`+`, Map(function(P, l) l * (P %*% beta), pen@mats, lam)))
}

#' @rdname penalty_gradient.AdditivePenalty
#' @name penalty_hessian.AdditivePenalty
#' @keywords internal
S7::method(penalty_hessian, AdditivePenalty) <- function(pen, beta, theta, ...) {
  lam <- unlist(theta[pen@params])
  Reduce(`+`, Map(function(P, l) l * P, pen@mats, lam))
}

#' @title Hyperparameter Derivatives of an Additive Penalty
#' @name penalty_grad_theta.AdditivePenalty
#'
#' @description
#' These methods return the three blocks of derivatives in the smoothing
#' parameters. `penalty_grad_theta()` returns one number per component,
#' `penalty_hess_theta()` one per unordered pair, and `penalty_cross()` one
#' coefficient vector per component.
#'
#' @details
#' Each parameter enters the value linearly through the quadratic form and
#' non-linearly through the log pseudo-determinant, so with \eqn{S^{+}} the
#' pseudo-inverse of the sum,
#'
#' \deqn{\frac{\partial\rho}{\partial\lambda_k}
#'     = \tfrac{1}{2}\beta^\top P_k\beta
#'       - \tfrac{1}{2}\operatorname{tr}(S^{+}P_k), \qquad
#'   \frac{\partial^2\rho}{\partial\lambda_k\partial\lambda_l}
#'     = \tfrac{1}{2}\operatorname{tr}(S^{+}P_kS^{+}P_l), \qquad
#'   \frac{\partial^2\rho}{\partial\beta\,\partial\lambda_k} = P_k\beta.}
#'
#' The second derivative has no term from the quadratic form, which is linear
#' in each parameter, and the mixed block has no term from the determinant.
#' `penalty_cross()` therefore takes no eigendecomposition. The first two
#' blocks read the traces from [additive_sum()], and their entries are `NaN`
#' at a setting that it reports as not resolvable.
#'
#' The Hessian is keyed by pairs as [penalty_hess_theta()] keys them: each key
#' joins two parameter names by an underscore (`lambda1_lambda2`), the
#' diagonal pairs come first and the upper off-diagonal pairs follow.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of `lambda1`, `lambda2`, ...
#' @param scale Read by the generic, which applies the chain rule onto the
#'   unconstrained scale after dispatch. These methods always return the
#'   parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` returns a list of one number per component,
#'   named `lambda1`, `lambda2`, ... `penalty_hess_theta()` returns a list of
#'   one number per unordered pair, keyed diagonals first. `penalty_cross()`
#'   returns a list of one numeric vector of length `pen@n_coef` per component.
#'
#' @examples
#' P1 <- crossprod(diff(diag(5)))
#' P2 <- crossprod(diff(diag(5), differences = 2))
#' pen <- additive_penalty(list(P1, P2))
#' th <- list(lambda1 = 2, lambda2 = 0.5)
#' b <- c(0.4, -1.1, 0.7, 0.2, -0.3)
#'
#' penalty_grad_theta(pen, b, th)
#' names(penalty_hess_theta(pen, b, th))
#'
#' # The mixed block is P_k beta, so it does not depend on the parameters.
#' identical(penalty_cross(pen, b, th),
#'           penalty_cross(pen, b, list(lambda1 = 99, lambda2 = 1e-3)))
#'
#' @seealso [penalty_value.AdditivePenalty()] for the quantity differentiated,
#'   [penalty_logpdet.AdditivePenalty()] for the determinant's own
#'   derivatives, [additive_sum()] for the traces that these read.
#' @keywords internal
S7::method(penalty_grad_theta, AdditivePenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    a <- additive_sum(pen, theta)
    stats::setNames(lapply(seq_along(pen@mats), function(k) {
      0.5 * sum(beta * as.numeric(pen@mats[[k]] %*% beta)) -
        0.5 * a$dlog[[k]]
    }), pen@params)
  }

#' @rdname penalty_grad_theta.AdditivePenalty
#' @name penalty_hess_theta.AdditivePenalty
#' @keywords internal
S7::method(penalty_hess_theta, AdditivePenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    a <- additive_sum(pen, theta)
    prs <- ptheta_pairs(pen@params)
    stats::setNames(lapply(names(prs), function(nm) {
      ij <- prs[[nm]]
      k <- match(ij[1], pen@params)
      l <- match(ij[2], pen@params)
      0.5 * a$d2log[k, l]
    }), names(prs))
  }

#' @rdname penalty_grad_theta.AdditivePenalty
#' @name penalty_cross.AdditivePenalty
#' @keywords internal
S7::method(penalty_cross, AdditivePenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    stats::setNames(lapply(pen@mats, function(P) as.numeric(P %*% beta)),
                    pen@params)
  }

#' @title Smoothness and Properness of an Additive Penalty
#' @name penalty_kinks.AdditivePenalty
#'
#' @description
#' `penalty_kinks()` returns `numeric(0)`, because the value is a quadratic
#' form and is smooth everywhere. `is_proper()` returns `TRUE` when the stored
#' rank equals the number of coefficients.
#'
#' @details
#' The value is a quadratic form, so [is_quadratic()] returns `TRUE` as well.
#' Smoothness does not imply a quadratic form: a Student t prior is smooth and
#' not quadratic.
#'
#' Properness compares the rank of the sum with the number of coefficients,
#' and the rank is fixed at construction. The null spaces of the two
#' components of an anisotropic tensor smooth have a nonzero intersection, so
#' such a penalty is improper: on a 4 by 4 grid penalized by second
#' differences along each margin, the rank is 12 out of 16.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param theta A named list of `lambda1`, `lambda2`, ... `penalty_kinks()`
#'   accepts it, and its result does not depend on it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` returns a numeric vector of length zero.
#'   `is_proper()` returns a single logical.
#'
#' @examples
#' P <- crossprod(diff(diag(4), differences = 2))
#' te <- additive_penalty(list(kronecker(diag(4), P), kronecker(P, diag(4))))
#'
#' penalty_kinks(te, list(lambda1 = 1, lambda2 = 1))
#' is_proper(te)
#' c(rank = penalty_rank(te), n_coef = te@n_coef)
#'
#' # A sum of full-rank components is proper.
#' is_proper(additive_penalty(list(diag(3), diag(c(2, 1, 1)))))
#'
#' @seealso [penalty_kinks()] and [is_proper()] for the generics,
#'   [is_quadratic.AdditivePenalty()] for the other predicate with a method
#'   on this class.
#' @keywords internal
S7::method(penalty_kinks, AdditivePenalty) <- function(pen, theta, ...) {
  numeric(0)
}

#' @rdname penalty_kinks.AdditivePenalty
#' @name is_proper.AdditivePenalty
#' @keywords internal
S7::method(is_proper, AdditivePenalty) <- function(pen, ...) {
  pen@p_rank == pen@n_coef
}

#' @title Marginal Quantities of an Additive Penalty
#' @name penalty_matrix.AdditivePenalty
#'
#' @description
#' These methods return the four quantities that a marginal criterion reads.
#' `penalty_matrix()` returns the weighted sum \eqn{S(\lambda) = \sum_k
#' \lambda_k P_k}, `penalty_rank()` the rank and `penalty_null_basis()` an
#' orthonormal basis of the null space shared by the components, both fixed at
#' construction, and `penalty_logpdet()` the log pseudo-determinant with its
#' first two derivatives.
#'
#' @details
#' With \eqn{S^{+}} the pseudo-inverse over the stored rank,
#'
#' \deqn{\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\,S
#'     = \operatorname{tr}(S^{+}P_k), \qquad
#'   \frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
#'     \log\mathrm{pdet}\,S = -\operatorname{tr}(S^{+}P_kS^{+}P_l).}
#'
#' Both formulas are exact. The value, the gradient and the Hessian of the log
#' pseudo-determinant come from the similarity transformation of
#' [additive_sum()], which takes the traces in the transformed coordinates and
#' not from a materialized \eqn{S^{+}}. All three are `NaN` at a setting that
#' [additive_sum()] reports as not resolvable. `penalty_matrix()` assembles
#' the sum directly, and `penalty_rank()` and `penalty_null_basis()` return
#' the stored values.
#'
#' `penalty_logpdet()` returns the same shape as
#' [penalty_logpdet.QuadraticPenalty()] and
#' [penalty_logpdet.StructuredPenalty()]: `grad` is a list keyed by
#' `pen@params` and `hess` a list keyed by the pairs that
#' [penalty_hess_theta()] uses, so code written for those branches reads
#' `grad$lambda1` here.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param theta A named list of `lambda1`, `lambda2`, ... Only
#'   `penalty_matrix()` and `penalty_logpdet()` read it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_matrix()` returns a symmetric base matrix of side
#'   `pen@n_coef`. `penalty_rank()` returns a single integer, and
#'   `penalty_null_basis()` an orthonormal basis of the null space shared by the
#'   components, a matrix with `pen@n_coef` rows and `pen@n_coef -
#'   penalty_rank(pen)` columns. `penalty_logpdet()` returns a list of `value`
#'   (a single number), `grad` (a list keyed by `pen@params`) and `hess` (a list
#'   keyed by the pairs that `penalty_hess_theta()` uses). Every entry of it is
#'   `NaN` at a setting that [additive_sum()] reports as not resolvable.
#'
#' @examples
#' P1 <- crossprod(diff(diag(5)))
#' P2 <- crossprod(diff(diag(5), differences = 2))
#' pen <- additive_penalty(list(P1, P2))
#' th <- list(lambda1 = 2, lambda2 = 0.5)
#'
#' penalty_rank(pen)
#' lp <- penalty_logpdet(pen, th)
#' str(lp)
#'
#' # The gradient is the trace of the pseudo-inverse against each component.
#' S <- penalty_matrix(pen, th)
#' e <- eigen(S, symmetric = TRUE)
#' k <- order(e$values, decreasing = TRUE)[seq_len(penalty_rank(pen))]
#' Sp <- e$vectors[, k, drop = FALSE] %*% (t(e$vectors[, k, drop = FALSE]) /
#'                                           e$values[k])
#' max(abs(unlist(lp$grad) - c(sum(Sp * P1), sum(Sp * P2))))
#'
#' # The null space of the sum is the intersection of the components', so it
#' # does not move with the hyperparameters.
#' nb <- penalty_null_basis(pen)
#' c(width = ncol(nb), annihilated = max(abs(S %*% nb)) < 1e-12)
#'
#' # And the value is the sum of the logarithms of the non-zero eigenvalues.
#' lp$value - sum(log(e$values[k]))
#'
#' @seealso [penalty_matrix()] and its siblings for the generics,
#'   [additive_sum()] for the decomposition these read,
#'   [penalty_matrix.QuadraticPenalty()] for the branch whose determinant is
#'   linear in one parameter.
#' @keywords internal
S7::method(penalty_matrix, AdditivePenalty) <- function(pen, theta, ...) {
  lam <- unlist(theta[pen@params])
  Reduce(`+`, Map(function(P, l) l * P, pen@mats, lam))
}

#' @rdname penalty_matrix.AdditivePenalty
#' @name penalty_rank.AdditivePenalty
#' @keywords internal
S7::method(penalty_rank, AdditivePenalty) <- function(pen, ...) pen@p_rank

#' @rdname penalty_matrix.AdditivePenalty
#' @name penalty_logpdet.AdditivePenalty
#' @keywords internal
S7::method(penalty_logpdet, AdditivePenalty) <- function(pen, theta, ...) {
  a <- additive_sum(pen, theta)
  prs <- ptheta_pairs(pen@params)
  list(
    value = a$logpdet,
    grad = stats::setNames(as.list(a$dlog), pen@params),
    hess = stats::setNames(lapply(names(prs), function(nm) {
      ij <- prs[[nm]]
      k <- match(ij[1], pen@params)
      l <- match(ij[2], pen@params)
      -a$d2log[k, l]
    }), names(prs))
  )
}
#' @rdname penalty_matrix.AdditivePenalty
#' @name penalty_null_basis.AdditivePenalty
#' @keywords internal
S7::method(penalty_null_basis, AdditivePenalty) <- function(pen, ...) {
  pen@null_basis
}


#' @title Quadratic Form of an Additive Penalty
#' @name is_quadratic.AdditivePenalty
#' @description
#' Returns `TRUE`. A sum of quadratic forms is a quadratic form, so this class
#' has methods for the matrix, the rank, the null basis and the log
#' pseudo-determinant that [is_quadratic()] gates, and a marginal criterion
#' can read them. The matrix moves with one hyperparameter per component,
#' whereas the plain quadratic branch has one scale, so the log
#' pseudo-determinant is not linear in the logarithm of any one parameter, and
#' [check_penalty()] compares its gradient with `numDeriv` where the plain
#' quadratic branch checks a slope.
#' @param pen An [AdditivePenalty()] object.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#' @return `TRUE`.
#' @seealso [penalty_matrix.AdditivePenalty()] for the quantities that this
#'   predicate gates.
#' @keywords internal
S7::method(is_quadratic, AdditivePenalty) <- function(pen, ...) TRUE


#' The Zero Matrix of an Additive Penalty's Side
#'
#' @description
#' Returns the zero matrix of side `pen@n_coef`: sparse when the penalty is
#' repeated blockwise, so that a third-order block added to it by a caller
#' is not allocated dense, and a base matrix otherwise.
#'
#' @param pen An [AdditivePenalty()] object.
#'
#' @return A square zero matrix.
#'
#' @keywords internal
additive_zero <- function(pen) {
  n <- as.integer(pen@n_coef)
  if (additive_blocks(pen) > 1L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(n, n)))
  }
  matrix(0, n, n)
}


#' @title Readable Quantities of an Additive Penalty
#' @name penalty_readable.AdditivePenalty
#'
#' @description
#' Returns `NULL` for a penalty built with `readable = "lambda"`, whose
#' smoothing parameters are the quantities reported, and the standard
#' deviations \eqn{\sigma_k = \lambda_k^{-1/2}} for one built with
#' `readable = "sd"`.
#'
#' @details
#' Each standard deviation depends only on its own parameter, so the Jacobian
#' is diagonal, with entries \eqn{-\tfrac{1}{2}\lambda_k^{-3/2}}, and each
#' interval is built on the log scale.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param theta A named list of `lambda1`, `lambda2`, ...
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `NULL`, or a list of `value` (named `sigma1`, `sigma2`, ...),
#'   `jacobian`, `transform` and `block`, as [penalty_readable()] documents.
#'
#' @examples
#' pen <- additive_penalty(list(diag(2), diag(c(1, 2))), readable = "sd")
#' penalty_readable(pen, list(lambda1 = 4, lambda2 = 0.25))
#'
#' @seealso [penalty_readable()] for the generic, [additive_penalty()] for
#'   the argument.
#' @keywords internal
S7::method(penalty_readable, AdditivePenalty) <- function(pen, theta, ...) {
  if (!identical(pen@readable, "sd")) return(NULL)
  lam <- unlist(theta[pen@params])
  K <- length(lam)
  nm <- paste0("sigma", seq_len(K))
  J <- diag(-0.5 * lam^(-1.5), K, K)
  dimnames(J) <- list(nm, pen@params)
  list(value = stats::setNames(lam^(-0.5), nm), jacobian = J,
       transform = rep("log", K), block = rep("sd", K))
}

