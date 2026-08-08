#' @include quadratic_penalty.R
NULL

#' @title S7 Class for a Sum of Quadratic Penalties
#'
#' @description
#' The class \code{\link{additive_penalty}} instantiates: several
#' quadratic penalties added together, each with a smoothing parameter of
#' its own.
#'
#' @inheritParams penalty
#' @param mats The list of penalty matrices, already carried through the
#'   map.
#' @param p_rank The rank of the sum, the same at every positive parameter
#'   value.
#'
#' @return An object of class \code{AdditivePenalty}.
#'
#' @seealso \code{\link{additive_penalty}}
#' @examples
#' S7::S7_inherits(additive_penalty(list(diag(3), diag(c(1, 0, 0)))),
#'                 AdditivePenalty)
#' @keywords internal
#' @export
AdditivePenalty <- S7::new_class(
  name = "AdditivePenalty",
  parent = penalty,
  properties = list(
    mats = S7::class_list,
    p_rank = S7::class_numeric
  )
)

#' @title Construct a Sum of Quadratic Penalties
#'
#' @description
#' The gaussian prior whose precision is a weighted sum of fixed matrices,
#' \deqn{\rho(\beta;\lambda) = \tfrac{1}{2}\beta^\top S(\lambda)\beta
#'   - \tfrac{1}{2}\log\mathrm{pdet}\,S(\lambda)
#'   + \tfrac{r}{2}\log 2\pi, \qquad
#'   S(\lambda) = \sum_{k} \lambda_k P_k,}
#' with one smoothing parameter per component. This is what a tensor-product
#' smooth needs: one parameter per margin, so that the fit may be rough in
#' one direction and smooth in another.
#'
#' @details
#' \code{\link{quadratic_penalty}} carries a single matrix and one scale,
#' which forces every direction to be smoothed alike. Here the components
#' keep their own parameters, and the quantities a marginal criterion reads
#' follow from one eigendecomposition of the sum:
#' \deqn{\frac{\partial}{\partial\lambda_k}\log\mathrm{pdet}\,S
#'   = \operatorname{tr}(S^{+}P_k), \qquad
#'   \frac{\partial^{2}}{\partial\lambda_k\partial\lambda_l}
#'   \log\mathrm{pdet}\,S = -\operatorname{tr}(S^{+}P_kS^{+}P_l).}
#'
#' \strong{The rank is not read off the sum.} Counting the eigenvalues of
#' \eqn{S(\lambda)} above a tolerance gives the right answer only while the
#' parameters are comparable: once they differ by orders of magnitude the
#' small contributions sink below any fixed tolerance and are counted as
#' zeros, so the rank appears to fall as the fit is smoothed. The null space
#' of a sum of positive semidefinite matrices is the intersection of their
#' null spaces, which does not depend on the parameters at all, so the rank
#' is fixed once at construction from the components stacked and
#' individually normalized.
#'
#' @param mats A list of symmetric positive semidefinite matrices of the
#'   same dimension.
#' @param map The matrix \eqn{D}, or \code{NULL} (default) for the
#'   identity.
#' @param link_lambda The link carrying each smoothing parameter to the
#'   unconstrained scale. Defaults to the log.
#' @param tol The relative tolerance below which an eigenvalue counts as
#'   zero.
#'
#' @return An object of class \code{\link{AdditivePenalty}}.
#'
#' @seealso \code{\link{quadratic_penalty}}
#'
#' @examples
#' # curvature in two directions, penalized separately
#' P1 <- crossprod(diff(diag(4), differences = 2))
#' pen <- additive_penalty(list(kronecker(diag(4), P1),
#'                              kronecker(P1, diag(4))))
#' pen@params
#' penalty_value(pen, rnorm(16), list(lambda1 = 1, lambda2 = 100))
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
    map <- as.matrix(map)
    if (nrow(map) != d) {
      stop(sprintf("'map' must have %d rows.", d), call. = FALSE)
    }
    mats <- lapply(mats, function(P) crossprod(map, P %*% map))
    d <- ncol(map)
  }

  # The rank is a property of the components, not of any one value of the
  # parameters: the null space of the sum is the intersection of the null
  # spaces, so the normalized components are added once and its rank read
  # there. Reading it off S(lambda) would make the rank fall as the
  # parameters spread apart.
  stacked <- Reduce(`+`, lapply(mats, function(P) P / max(abs(P))))
  ev <- eigen(stacked, symmetric = TRUE, only.values = TRUE)$values
  r <- sum(ev > tol * max(ev))

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
    p_rank = r
  )
}

#' The Weighted Sum and Its Pseudo-Inverse
#'
#' @description
#' \eqn{S(\lambda)} and the pseudo-inverse the log pseudo-determinant and
#' its derivatives are written in, from one eigendecomposition.
#'
#' @param pen An \code{AdditivePenalty} object.
#' @param theta The aligned hyperparameter list.
#'
#' @return A list with the matrix \code{S}, its pseudo-inverse \code{Sp},
#'   and the log pseudo-determinant \code{logpdet}.
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

#' @title Additive Penalty Methods
#' @name penalty_value.AdditivePenalty
#' @description The closed forms; see \code{\link{additive_penalty}}.
#' @param pen An \code{AdditivePenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of smoothing parameters.
#' @param scale Handled by the generic.
#' @param ... Unused.
#' @return See the generic pages.
#' @keywords internal
S7::method(penalty_value, AdditivePenalty) <- function(pen, beta, theta, ...) {
  a <- additive_sum(pen, theta)
  0.5 * sum(beta * (a$S %*% beta)) - 0.5 * a$logpdet +
    pen@p_rank / 2 * log(2 * pi)
}

#' @rdname penalty_value.AdditivePenalty
#' @name penalty_gradient.AdditivePenalty
#' @keywords internal
S7::method(penalty_gradient, AdditivePenalty) <- function(pen, beta, theta, ...) {
  lam <- unlist(theta[pen@params])
  as.numeric(Reduce(`+`, Map(function(P, l) l * (P %*% beta), pen@mats, lam)))
}

#' @rdname penalty_value.AdditivePenalty
#' @name penalty_hessian.AdditivePenalty
#' @keywords internal
S7::method(penalty_hessian, AdditivePenalty) <- function(pen, beta, theta, ...) {
  lam <- unlist(theta[pen@params])
  Reduce(`+`, Map(function(P, l) l * P, pen@mats, lam))
}

#' @rdname penalty_value.AdditivePenalty
#' @name penalty_grad_theta.AdditivePenalty
#' @keywords internal
S7::method(penalty_grad_theta, AdditivePenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    a <- additive_sum(pen, theta)
    stats::setNames(lapply(pen@mats, function(P) {
      0.5 * sum(beta * (P %*% beta)) - 0.5 * sum(a$Sp * P)
    }), pen@params)
  }

#' @rdname penalty_value.AdditivePenalty
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

#' @rdname penalty_value.AdditivePenalty
#' @name penalty_cross.AdditivePenalty
#' @keywords internal
S7::method(penalty_cross, AdditivePenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    stats::setNames(lapply(pen@mats, function(P) as.numeric(P %*% beta)),
                    pen@params)
  }

#' @rdname penalty_value.AdditivePenalty
#' @name penalty_kinks.AdditivePenalty
#' @keywords internal
S7::method(penalty_kinks, AdditivePenalty) <- function(pen, theta, ...) {
  numeric(0)
}

#' @rdname penalty_value.AdditivePenalty
#' @name is_proper.AdditivePenalty
#' @keywords internal
S7::method(is_proper, AdditivePenalty) <- function(pen, ...) {
  pen@p_rank == pen@n_coef
}

#' @rdname penalty_value.AdditivePenalty
#' @name penalty_matrix.AdditivePenalty
#' @keywords internal
S7::method(penalty_matrix, AdditivePenalty) <- function(pen, theta, ...) {
  lam <- unlist(theta[pen@params])
  Reduce(`+`, Map(function(P, l) l * P, pen@mats, lam))
}

#' @rdname penalty_value.AdditivePenalty
#' @name penalty_rank.AdditivePenalty
#' @keywords internal
S7::method(penalty_rank, AdditivePenalty) <- function(pen, ...) pen@p_rank

#' @rdname penalty_value.AdditivePenalty
#' @name penalty_logpdet.AdditivePenalty
#' @keywords internal
S7::method(penalty_logpdet, AdditivePenalty) <- function(pen, theta, ...) {
  a <- additive_sum(pen, theta)
  list(value = a$logpdet,
       grad = vapply(pen@mats, function(P) sum(a$Sp * P), numeric(1)),
       hess = outer(seq_along(pen@mats), seq_along(pen@mats),
                    Vectorize(function(k, l)
                      -sum(t(a$Sp %*% pen@mats[[k]]) * (a$Sp %*% pen@mats[[l]])))))
}
