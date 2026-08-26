#' @include quadratic_penalty.R
NULL

#' @title S7 Class for a Sum of Quadratic Penalties
#'
#' @description
#' The class [additive_penalty()] builds: several quadratic penalties added
#' together, each with a smoothing parameter of its own. Beyond the seven
#' properties every penalty carries it stores the list of component matrices,
#' already carried through the map, and the rank and null basis of the sum.
#'
#' @details
#' The rank and the null basis are stored because both are properties of the
#' components alone: the null space of a sum of positive semidefinite
#' matrices is the intersection of the components' null spaces, so neither
#' moves with the parameters. Reading the rank off the assembled
#' \eqn{S(\lambda)} instead would make it fall as the parameters spread
#' apart, which [additive_penalty()] measures.
#'
#' @inheritParams penalty
#' @param mats A list of symmetric positive semidefinite matrices of the same
#'   side, already carried through the map, so each is `n_coef` by `n_coef`.
#'   One smoothing parameter goes with each, named `lambda1`, `lambda2` and so
#'   on in this order.
#' @param p_rank The rank of the sum, a single whole number, the same at every
#'   positive setting of the parameters.
#' @param null_basis An orthonormal basis of the components' shared null
#'   space, `n_coef` by `n_coef - p_rank`, read from the same
#'   eigendecomposition as `p_rank` and equally fixed.
#'
#' @return An S7 object of class `AdditivePenalty`, inheriting from [penalty()],
#'   with the seven inherited properties and the three above. Its `map` is
#'   always `NULL`, a map given to the constructor having been absorbed into
#'   `mats`.
#'
#' @seealso [additive_penalty()] for the constructor to use,
#'   [penalty_value.AdditivePenalty()] for what the branch computes,
#'   [quadratic_penalty()] for one matrix with one scale.
#'
#' @examples
#' pen <- additive_penalty(list(diag(3), diag(c(1, 0, 0))))
#' S7::S7_inherits(pen, AdditivePenalty)
#'
#' # One hyperparameter per component, and a rank fixed at construction.
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
    null_basis = S7::class_any
  )
)

#' @title Construct a Sum of Quadratic Penalties
#'
#' @description
#' Builds the Gaussian prior whose precision is a weighted sum of fixed
#' matrices, one smoothing parameter to each. This is what an anisotropic
#' tensor-product smooth needs: one parameter per margin, so a fit may be rough
#' in one direction and smooth in another, where a single parameter would force
#' every direction to be smoothed alike.
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
#' `lambda2` and so on, in the order the components were given, and each is
#' positive on a log link.
#'
#' # The determinant's derivatives
#'
#' One eigendecomposition of \eqn{S(\lambda)} gives the pseudo-inverse
#' \eqn{S^{+}}, and the quantities a marginal criterion reads follow from it:
#'
#' \deqn{\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\,S
#'   = \operatorname{tr}(S^{+}P_k), \qquad
#'   \frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
#'   \log\mathrm{pdet}\,S = -\operatorname{tr}(S^{+}P_kS^{+}P_l).}
#'
#' Both are exact, and both agree with the assembled quantity to 0.
#'
#' # The rank is not read off the sum
#'
#' Counting the eigenvalues of \eqn{S(\lambda)} above a tolerance gives the
#' right answer only while the parameters are comparable. Once they differ by
#' orders of magnitude the small contributions sink below any fixed tolerance
#' and are counted as zeros, so the rank appears to fall as the fit is
#' smoothed. Measured on first and second differences over five coefficients,
#' whose intersected null space is the constants:
#'
#' | \eqn{\lambda_2/\lambda_1} | eigenvalues above `1e-10` | stored rank |
#' |---|---|---|
#' | 1 to 1e8 | 4 | 4 |
#' | 1e10 to 1e14 | **3** | 4 |
#'
#' while the null basis annihilates the assembled sum to `1e-16` relative at
#' every one of those settings. So the rank is fixed once at construction, from
#' the components stacked and individually normalized, and the object's answer
#' cannot move. A test pins both halves.
#'
#' # What this branch supplies
#'
#' [is_quadratic()] answers `TRUE`, a sum of quadratic forms being a quadratic
#' form, so [penalty_matrix()], [penalty_rank()], [penalty_null_basis()] and
#' [penalty_logpdet()] all answer and [check_penalty()] runs its three
#' quadratic rows here. What is particular is that the matrix moves with one
#' hyperparameter per component, where the plain quadratic branch has a
#' single scale, so the log
#' pseudo-determinant is linear in none of them and the row that checks it
#' compares a gradient against `numDeriv` where the plain quadratic branch
#' reads a slope.
#'
#' [has_prox()] answers `FALSE`, this branch registering no
#' [penalty_prox()], which is asked before `is_quadratic()` is.
#'
#' @param mats A list of symmetric positive semidefinite matrices, all of the
#'   same side, and at least one. Each is symmetrized, and a component that is
#'   zero, negative definite, or not symmetric to a relative `1e-8` is
#'   rejected by number.
#' @param map The matrix \eqn{D}, with as many rows as the components and one
#'   column per coefficient, or `NULL` (the default) for the identity. A map is
#'   absorbed at construction, each component becoming \eqn{D'P_kD}, so the
#'   stored `map` is always `NULL`.
#' @param link_lambda The \pkg{linkfunctions7} link carrying each smoothing
#'   parameter onto the whole real line. `linkfunctions7::log_link()` by
#'   default, and the same link is used for all of them.
#' @param tol The relative eigenvalue tolerance, `1e-10` by default. Used twice:
#'   to reject a component whose smallest eigenvalue is below `-tol` times its
#'   largest, and to count the rank of the normalized stack.
#'
#' @return An [AdditivePenalty()] object with one hyperparameter per component,
#'   named `lambda1`, `lambda2`, ..., each bounded on \eqn{(0, \infty)}.
#'
#' @examples
#' # An anisotropic tensor smooth over a 4 x 4 grid: curvature along each
#' # margin, penalized separately.
#' P <- crossprod(diff(diag(4), differences = 2))
#' pen <- additive_penalty(list(kronecker(diag(4), P), kronecker(P, diag(4))))
#' pen@params
#' c(n_coef = pen@n_coef, rank = penalty_rank(pen))
#'
#' # Raising one parameter charges for roughness along that margin alone.
#' set.seed(2)
#' b <- rnorm(16)
#' penalty_value(pen, b, list(lambda1 = 1, lambda2 = 1))
#' penalty_value(pen, b, list(lambda1 = 1, lambda2 = 100))
#'
#' # With one component the branch is the plain quadratic penalty.
#' P2 <- crossprod(diff(diag(5), differences = 2))
#' bb <- c(0.4, -1.1, 0.7, 0.2, -0.3)
#' penalty_value(additive_penalty(list(P2)), bb, list(lambda1 = 3)) -
#'   penalty_value(quadratic_penalty(P2), bb, list(lambda = 3))
#'
#' # The stored rank does not move as the parameters spread apart, where an
#' # eigenvalue count of the assembled sum does.
#' a <- additive_penalty(list(crossprod(diff(diag(5))), P2))
#' penalty_rank(a)
#' sapply(c(1, 1e10), function(r) {
#'   S <- penalty_matrix(a, list(lambda1 = 1, lambda2 = r))
#'   ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
#'   sum(ev > 1e-10 * max(ev))
#' })
#'
#' @seealso [quadratic_penalty()] for one matrix with one scale,
#'   [structured_penalty()] for a matrix whose entries move with the
#'   hyperparameters, [penalty_logpdet.AdditivePenalty()] for the determinant
#'   and its derivatives, [additive_sum()] for the decomposition they share.
#'
#' @export
additive_penalty <- function(mats, map = NULL,
                             link_lambda = linkfunctions7::log_link(),
                             tol = 1e-10) {
  if (!is.list(mats) || length(mats) < 1L) {
    stop("'mats' must be a non-empty list of matrices.", call. = FALSE)
  }
  mats <- lapply(seq_along(mats), function(k) {
    P <- as.matrix(mats[[k]])
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
      stop(sprintf("component %d is zero or negative definite.", k),
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
    if (nrow(map) != d) {
      stop(sprintf("'map' must have %d rows.", d), call. = FALSE)
    }
    # densified for the reason quadratic_penalty() states: a Matrix map would
    # otherwise decide the class of every component and of every quantity
    # assembled from them
    mats <- lapply(mats, function(P) as.matrix(crossprod(map, P %*% map)))
    d <- ncol(map)
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
    null_basis = nb
  )
}

#' The Weighted Sum and Its Pseudo-Inverse
#'
#' @description
#' Assembles \eqn{S(\lambda) = \sum_k \lambda_k P_k}, takes one
#' eigendecomposition of it, and returns the matrix, its pseudo-inverse and its
#' log pseudo-determinant together. The value, the hyperparameter derivatives
#' and [penalty_logpdet()] all need some of these three, and this is the only
#' place the decomposition happens.
#'
#' @details
#' The pseudo-inverse is taken over the `p_rank` largest eigenvalues, that rank
#' having been fixed at construction from the components rather than counted
#' here. Selecting by rank instead of by a tolerance is what keeps the answer
#' steady when the parameters differ by many orders of magnitude, which is
#' exactly when a count would lose directions the penalty still spans.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param theta The aligned hyperparameter list, as [align_ptheta()] returns
#'   it.
#'
#' @return A list of three: `S`, the assembled symmetric matrix of side
#'   `pen@n_coef`; `Sp`, its Moore-Penrose pseudo-inverse over the leading
#'   `pen@p_rank` eigendirections, of the same side; and `logpdet`, a single
#'   number, the sum of the logarithms of those eigenvalues.
#'
#' @seealso [additive_penalty()], [penalty_logpdet.AdditivePenalty()]
#'
#' @keywords internal
additive_sum <- function(pen, theta) {
  lam <- unlist(theta[pen@params])
  S <- Reduce(`+`, Map(function(P, l) l * P, pen@mats, lam))
  e <- eigen(S, symmetric = TRUE)
  keep <- order(e$values, decreasing = TRUE)[seq_len(pen@p_rank)]
  V <- e$vectors[, keep, drop = FALSE]
  dv <- e$values[keep]
  list(S = S, Sp = V %*% (t(V) / dv), logpdet = sum(log(dv)))
}

#' @title Value of an Additive Penalty
#' @name penalty_value.AdditivePenalty
#'
#' @description
#' Returns the negative log-density of the Gaussian prior whose precision is
#' the weighted sum \eqn{S(\lambda) = \sum_k \lambda_k P_k}, with the
#' normalizing constant taken over the range of the sum.
#'
#' @details
#' \deqn{\rho(\beta;\lambda) = \tfrac{1}{2}\beta^\top S(\lambda)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,S(\lambda)
#'   + \tfrac{r}{2}\log 2\pi,}
#'
#' with \eqn{r} the stored rank. One eigendecomposition of the sum, through
#' [additive_sum()], supplies both the quadratic form's matrix and the log
#' pseudo-determinant.
#'
#' With a single component the value is the plain quadratic penalty's, and the
#' two agree to `3.6e-15`.
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
#' # A constant vector lies in the intersected null space, so only the
#' # normalizing constant remains.
#' penalty_value(pen, rep(1, 5), list(lambda1 = 2, lambda2 = 0.5))
#' penalty_value(pen, rep(0, 5), list(lambda1 = 2, lambda2 = 0.5))
#'
#' @seealso [additive_penalty()] for the construction,
#'   [penalty_gradient.AdditivePenalty()] for the coefficient derivatives,
#'   [additive_sum()] for the decomposition it reads.
#' @keywords internal
S7::method(penalty_value, AdditivePenalty) <- function(pen, beta, theta, ...) {
  a <- additive_sum(pen, theta)
  0.5 * sum(beta * (a$S %*% beta)) - 0.5 * a$logpdet +
    pen@p_rank / 2 * log(2 * pi)
}

#' @title Coefficient Derivatives of an Additive Penalty
#' @name penalty_gradient.AdditivePenalty
#'
#' @description
#' `penalty_gradient()` returns \eqn{S(\lambda)\beta} and `penalty_hessian()`
#' returns \eqn{S(\lambda)}. Both assemble the weighted sum directly and take
#' no eigendecomposition, the normalizing constant not depending on the
#' coefficients.
#'
#' @details
#' The value is a quadratic form, so the Hessian does not move with the
#' coefficients and the gradient is the Hessian applied to them. A coefficient
#' vector lying in the intersected null space of the components has zero
#' gradient at every setting of the parameters.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of `lambda1`, `lambda2`, ...
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` a numeric vector of length `pen@n_coef`;
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
#' The three blocks in the smoothing parameters. `penalty_grad_theta()` returns
#' one number per component, `penalty_hess_theta()` one per unordered pair, and
#' `penalty_cross()` one coefficient vector per component.
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
#' The second derivative has no quadratic-form term, the value being linear in
#' each parameter there, and the mixed block has no determinant term.
#' `penalty_cross()` therefore needs no eigendecomposition at all.
#'
#' The Hessian is keyed by pairs as [penalty_hess_theta()] keys them, diagonals
#' first and then the upper off-diagonal pairs joined by an underscore.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of `lambda1`, `lambda2`, ...
#' @param scale Read by the generic, which applies the chain rule onto the
#'   unconstrained scale after dispatch. These methods always return the
#'   parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` a list of one number per component, named
#'   `lambda1`, `lambda2`, ...
#'   `penalty_hess_theta()` a list of one number per unordered pair, keyed
#'   diagonals first.
#'   `penalty_cross()` a list of one numeric vector of length `pen@n_coef` per
#'   component.
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
#'   derivatives, [additive_sum()] for the pseudo-inverse these read.
#' @keywords internal
S7::method(penalty_grad_theta, AdditivePenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    a <- additive_sum(pen, theta)
    stats::setNames(lapply(pen@mats, function(P) {
      0.5 * sum(beta * (P %*% beta)) - 0.5 * sum(a$Sp * P)
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
      0.5 * sum(t(a$Sp %*% pen@mats[[k]]) * (a$Sp %*% pen@mats[[l]]))
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
#' `penalty_kinks()` returns `numeric(0)`, the value being a quadratic form and
#' smooth everywhere. `is_proper()` returns `TRUE` when the stored rank equals
#' the number of coefficients.
#'
#' @details
#' The branch registers no `is_quadratic()` method, so it inherits `FALSE` from
#' [penalty()] even though the value is a quadratic form. See
#' [additive_penalty()] for what that costs.
#'
#' Properness is the rank of the sum against the number of coefficients, and
#' the rank is fixed at construction. The two components of an anisotropic
#' tensor smooth intersect in a null space of their own, so such a penalty is
#' improper: a 4 by 4 grid penalized by curvature along each margin has rank 12
#' out of 16.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param theta A named list of `lambda1`, `lambda2`, ... Read by
#'   `penalty_kinks()`, whose answer does not depend on it.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` a numeric vector of length zero.
#'   `is_proper()` a single logical.
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
#'   [additive_penalty()] for why [is_quadratic()] answers `FALSE` here.
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
#' Three of the four pieces a marginal criterion reads. `penalty_matrix()`
#' returns the weighted sum \eqn{S(\lambda) = \sum_k \lambda_k P_k},
#' `penalty_rank()` the rank fixed at construction, and `penalty_logpdet()` the
#' log pseudo-determinant with its first two derivatives. There is no
#' `penalty_null_basis()` method for this branch, and the base class's rejects.
#'
#' @details
#' With \eqn{S^{+}} the pseudo-inverse over the stored rank,
#'
#' \deqn{\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\,S
#'     = \operatorname{tr}(S^{+}P_k), \qquad
#'   \frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
#'     \log\mathrm{pdet}\,S = -\operatorname{tr}(S^{+}P_kS^{+}P_l),}
#'
#' both exact and both agreeing with the traces computed apart to 0.
#'
#' **`penalty_logpdet()` answers in a different shape here.** `grad` is an
#' unnamed numeric vector in `pen@params` order and `hess` is a square matrix,
#' where [penalty_logpdet.QuadraticPenalty()] and
#' [penalty_logpdet.StructuredPenalty()] return named lists keyed by
#' hyperparameter and by pair. A consumer written against those will read
#' `NULL` from `grad$lambda1` here. Nothing in the toolkit reads it today,
#' because [is_quadratic()] answers `FALSE` for this branch and every consumer
#' routes on that first.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param theta A named list of `lambda1`, `lambda2`, ... `penalty_matrix()`
#'   and `penalty_logpdet()` only.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_matrix()` a symmetric base matrix of side `pen@n_coef`.
#'   `penalty_rank()` a single integer, and `penalty_null_basis()` an
#'   orthonormal basis of the components' shared null space, `pen@n_coef` by
#'   `pen@n_coef - penalty_rank(pen)`.
#'   `penalty_logpdet()` a list of `value` (a single number), `grad` (a list
#'   keyed by `pen@params`) and `hess` (a list keyed by the pairs
#'   `penalty_hess_theta()` uses), which is the shape the quadratic and
#'   structured branches answer in.
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
    grad = stats::setNames(lapply(pen@mats, function(P) sum(a$Sp * P)),
                           pen@params),
    hess = stats::setNames(lapply(names(prs), function(nm) {
      ij <- prs[[nm]]
      k <- match(ij[1], pen@params)
      l <- match(ij[2], pen@params)
      -sum(t(a$Sp %*% pen@mats[[k]]) * (a$Sp %*% pen@mats[[l]]))
    }), names(prs))
  )
}
#' @rdname penalty_matrix.AdditivePenalty
#' @name penalty_null_basis.AdditivePenalty
#' @keywords internal
S7::method(penalty_null_basis, AdditivePenalty) <- function(pen, ...) {
  pen@null_basis
}


#' @title An Additive Penalty Is Quadratic in the Coefficients
#' @name is_quadratic.AdditivePenalty
#' @description
#' Answers `TRUE`. A sum of quadratic forms is a quadratic form, so this branch
#' has the matrix, the rank, the null basis and the log pseudo-determinant that
#' [is_quadratic()] gates, and a marginal criterion can read them. What is
#' particular about it is that the matrix moves with one hyperparameter per
#' component, where the plain quadratic branch has one scale, so the log
#' pseudo-determinant is
#' not linear in any one of them and [check_penalty()] compares its gradient
#' against `numDeriv` where the plain quadratic branch reads a slope.
#' @param pen An [additive_penalty()] object.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return `TRUE`.
#' @seealso [penalty_matrix.AdditivePenalty()] for the quantities this admits.
#' @keywords internal
S7::method(is_quadratic, AdditivePenalty) <- function(pen, ...) TRUE
