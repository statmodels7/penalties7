#' @include generics.R quadratic_penalty.R additive_penalty.R structured_penalty.R distrib_penalty.R scad_mcp.R
NULL

# What a marginal criterion asks of a penalty, beyond the value, the gradient,
# the Hessian and the mixed block.
#
# Estimating the hyperparameters by a marginal likelihood needs the derivative
# of the LAPLACE approximation, whose determinant is |H + S| with S the
# penalty's Hessian in the coefficients. Differentiating it in the
# hyperparameters therefore needs dS/dtheta, and differentiating a second time
# needs d2S/dtheta2 and d3rho/dbeta dtheta2. None of those is expressible from
# the value, the gradient and the Hessian: they are new quantities, and they
# belong here rather than in whatever consumes them, so that a penalty written
# later answers them by writing three methods and works in a marginal
# criterion with no change anywhere else.
#
# The base class rejects. A penalty that cannot answer says so, rather than
# being read through a test of its behavior -- a consumer that measured
# whether a Hessian happened to be linear in the hyperparameters would be
# guessing at a property the penalty knows.

#' The Derivative of the Coefficient Hessian in the Hyperparameters
#'
#' @description
#' Returns \eqn{\partial S/\partial\theta_m = \partial^3\rho/\partial\beta^2\,
#' \partial\theta_m}, one matrix per hyperparameter, which describes how the
#' curvature of the penalty in the coefficients changes with each
#' hyperparameter. A marginal likelihood criterion needs this third-order
#' quantity, and the value, the gradient and the Hessian do not determine it.
#'
#' @details
#' # Role in a marginal criterion
#'
#' Estimating hyperparameters by a marginal likelihood means differentiating a
#' Laplace approximation, whose determinant is \eqn{\lvert H + S\rvert} with
#' \eqn{S} the Hessian of the penalty in the coefficients. The derivative of
#' \eqn{\log\lvert H+S\rvert} in \eqn{\theta_m} contains
#' \eqn{\partial S/\partial\theta_m}, which is a third derivative of
#' \eqn{\rho}. The derivatives are taken in \eqn{\theta} on its own scale, and
#' no link is applied.
#'
#' # The value on each branch
#'
#' | branch | \eqn{\partial S/\partial\theta_m} |
#' |---|---|
#' | [quadratic_penalty()] | \eqn{D'PD}, free of both \eqn{\lambda} and \eqn{\beta} |
#' | [additive_penalty()] | the component \eqn{D'P_kD} of each smoothing parameter, which is \eqn{P_k} under the identity map |
#' | [structured_penalty()] | the derivative `param_d1` of the matrix parameter |
#' | [distrib_penalty()], univariate parent | \eqn{-D'\mathrm{diag}(\partial^3\ell/\partial y^2\partial\theta_m)D}, from [distributions7::distrib_cross2_y()] |
#' | [distrib_penalty()], multivariate parent | the same with a block-diagonal middle matrix |
#' | [scad_penalty()], [mcp_penalty()] | an error, from the base class |
#'
#' The quadratic and additive rows are exact. The structured row is exact for
#' the structures that \pkg{parameters7} ships, and the separable rows are
#' exact when the parent supplies [distributions7::distrib_cross2_y()] in
#' closed form; otherwise the numerical fallbacks of those packages apply. On
#' the base class the method signals an error that names the penalty, so a
#' penalty class written later is usable in a marginal criterion only if it
#' supplies `penalty_dhessian()`, [penalty_d2hessian()] and
#' [penalty_dcross()].
#'
#' A separable penalty with a **kink** signals an error as well, whatever its
#' parent. The third derivative does not exist at the kink, and the mode at
#' which a marginal criterion is expanded is where a kinked penalty places
#' coefficients.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`; the generic signals an
#'   error for any other length. Only the separable branch reads its values,
#'   because the derivatives of the parent depend on them; the quadratic,
#'   additive and structured branches ignore them.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same. The generic checks that every hyperparameter is
#'   present as a single number inside its bounds.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A named list of one matrix per hyperparameter, each of side
#'   `pen@n_coef`, named by `pen@params`. The quadratic branch returns the
#'   stored matrix in the storage that the penalty keeps, so a penalty built
#'   with `blocks > 1` gives a `dgCMatrix`; the other branches return base
#'   matrices.
#'
#' @seealso [penalty_hessian()] for the quantity differentiated,
#'   [penalty_d2hessian()] and [penalty_dcross()] for the other two quantities
#'   that a marginal criterion requires, [beta_quadratic()] for the predicate
#'   that indicates whether a third \eqn{\beta}-derivative is needed.
#'
#' @examples
#' # A quadratic penalty's Hessian is lambda D'PD, so the derivative is D'PD
#' # and does not change with lambda.
#' pen <- quadratic_penalty(crossprod(diff(diag(3))))
#' b <- c(1, -0.5, 0.3)
#' penalty_dhessian(pen, b, list(lambda = 2))$lambda
#' identical(penalty_dhessian(pen, b, list(lambda = 2)),
#'           penalty_dhessian(pen, b, list(lambda = 99)))
#'
#' # An additive penalty returns its components, one per parameter.
#' add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
#' names(penalty_dhessian(add, b, list(lambda1 = 2, lambda2 = 0.5)))
#'
#' # A penalty with a kink has no such quantity, and the call is an error.
#' try(penalty_dhessian(lasso_penalty(n_coef = 3), b, list(lambda = 1)))
#'
#' @export
penalty_dhessian <- S7::new_generic("penalty_dhessian", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })


#' The Second Derivative of the Coefficient Hessian in the Hyperparameters
#'
#' @description
#' Returns \eqn{\partial^2 S/\partial\theta_m\partial\theta_l =
#' \partial^4\rho/\partial\beta^2\partial\theta_m\partial\theta_l}, one matrix
#' per unordered pair. The second derivative of a marginal criterion in the
#' hyperparameters contains it, and so does an exact outer Hessian.
#'
#' @details
#' It is **exactly zero** for every penalty whose Hessian is linear in its
#' hyperparameters, which is the quadratic branch, where
#' \eqn{S = \lambda D'PD}, and the additive branch, where
#' \eqn{S = \sum_k \lambda_k D'P_kD}. The structured branch reads the second
#' derivative `param_d2` of the matrix parameter, and the separable branch the
#' second \eqn{\theta}-derivative of the response curvature of the parent,
#' through [distributions7::distrib_hess_y_hess()], carried through the map as
#' in [penalty_dhessian()]. The derivatives are taken in \eqn{\theta} on its
#' own scale.
#'
#' The keys are those of [penalty_hess_theta()]: the diagonal pairs first in
#' `pen@params` order, then the upper off-diagonal pairs, each pair joined by
#' an underscore. Each unordered pair has one key, written in `pen@params`
#' order (`sigma_nu` and not `nu_sigma` when `sigma` comes first), and the
#' reversed name is not an element of the list.
#'
#' A penalty with a kink signals an error, and so does the base class.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`; the generic signals an
#'   error for any other length. Only the separable branch reads its values.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, checked by the generic as in
#'   [penalty_dhessian()].
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A named list of one matrix per unordered pair of hyperparameters,
#'   each of side `pen@n_coef`, keyed diagonals first. A penalty with \eqn{p}
#'   hyperparameters gives \eqn{p(p+1)/2} entries.
#'
#' @seealso [penalty_dhessian()] for the first order, [penalty_dcross()] for
#'   the mixed one, [penalty_hess_theta()] for the keying.
#'
#' @examples
#' b <- c(1, -0.5, 0.3)
#'
#' # Zero for a quadratic penalty, whose Hessian is linear in lambda.
#' pen <- quadratic_penalty(diag(3))
#' penalty_d2hessian(pen, b, list(lambda = 2))
#'
#' # Zero for an additive one as well, keyed diagonals first over three pairs.
#' add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
#' names(penalty_d2hessian(add, b, list(lambda1 = 2, lambda2 = 0.5)))
#'
#' # Not zero for a structured penalty, whose matrix is not linear in its
#' # free values.
#' s <- structured_penalty(parameters7::log_cholesky(2))
#' th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
#' penalty_d2hessian(s, c(1, -0.5), th)$log_L1_log_L1
#'
#' @export
penalty_d2hessian <- S7::new_generic("penalty_d2hessian", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })


#' The Derivative of the Mixed Block in the Hyperparameters
#'
#' @description
#' Returns \eqn{\partial^3\rho/\partial\beta\,\partial\theta_m\partial\theta_l},
#' one coefficient vector per unordered pair. [penalty_cross()] gives the
#' derivative of the coefficient gradient in one hyperparameter, and this
#' generic gives the derivative of that quantity in a second hyperparameter.
#'
#' @details
#' A marginal criterion is evaluated at the penalized mode, and the mode moves
#' with the hyperparameters. The second derivative of the criterion therefore
#' contains the second derivative of the mode, and this quantity is the
#' contribution of the penalty to it.
#'
#' It is **exactly zero** wherever the penalty is quadratic in the coefficients
#' with a Hessian linear in the hyperparameters, which covers the quadratic and
#' additive branches. The structured branch returns
#' \eqn{(\partial^2\Omega/\partial\theta_m\partial\theta_l)\beta}, computed
#' from [penalty_d2hessian()], and the separable branch returns
#' \eqn{-D'(\partial^3\ell/\partial y\,\partial\theta_m\partial\theta_l)}, from
#' [distributions7::distrib_grad_y_hess()] of the parent evaluated at
#' \eqn{D\beta}.
#'
#' The keys are those of [penalty_hess_theta()], diagonal pairs first.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`; the generic signals an
#'   error for any other length. The structured and separable branches read
#'   its values; the quadratic and additive branches return zero for any
#'   value.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, checked by the generic as in
#'   [penalty_dhessian()].
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A named list of one numeric vector of length `pen@n_coef` per
#'   unordered pair, keyed diagonals first.
#'
#' @seealso [penalty_cross()] for the first order,
#'   [penalty_dhessian()] and [penalty_d2hessian()] for the other two
#'   quantities that a marginal criterion requires.
#'
#' @examples
#' b <- c(1, -0.5, 0.3)
#'
#' # Zero on the quadratic branch.
#' penalty_dcross(quadratic_penalty(diag(3)), b, list(lambda = 2))
#'
#' # And not on the structured one, where it is d2 Omega / dtheta^2 times
#' # the coefficients.
#' s <- structured_penalty(parameters7::log_cholesky(2))
#' th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
#' bb <- c(1, -0.5)
#' penalty_dcross(s, bb, th)$log_L1_log_L1
#' drop(penalty_d2hessian(s, bb, th)$log_L1_log_L1 %*% bb)
#'
#' @export
penalty_dcross <- S7::new_generic("penalty_dcross", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })


#' Is a Penalty Quadratic in the Coefficients?
#'
#' @description
#' `TRUE` when \eqn{\partial^3\rho/\partial\beta^3} is exactly zero, so that a
#' consumer can skip that third derivative. The value comes from a method of
#' the class of the penalty, so a consumer does not have to test the property
#' numerically.
#'
#' @details
#' # Three independent properties
#'
#' This predicate differs from [is_quadratic()], which reports whether the
#' penalty is a quadratic form whose matrix, rank and log pseudo-determinant
#' are available through [penalty_matrix()] and its siblings. The two
#' predicates and the linearity of \eqn{S} in \eqn{\theta} are independent
#' properties:
#'
#' | penalty | `is_quadratic()` | `beta_quadratic()` | \eqn{S} linear in \eqn{\theta} |
#' |---|---|---|---|
#' | [quadratic_penalty()] | `TRUE` | `TRUE` | yes |
#' | [structured_penalty()] | `TRUE` | `TRUE` | no |
#' | [additive_penalty()] | `TRUE` | `TRUE` | yes |
#' | a Gaussian [distrib_penalty()] | `FALSE` | `TRUE` | no |
#' | [heavy_penalty()] | `FALSE` | `FALSE` | no |
#'
#' # The method of the separable branch
#'
#' The log-density is quadratic in the response exactly when its second
#' derivative in the response does not depend on the response. The method
#' therefore evaluates [distributions7::distrib_hess_y()] of the parent at four
#' points and tests whether the values agree across the points with
#' [all.equal()] at a tolerance of `1e-12`, which is a relative tolerance. No
#' third derivative is computed. The page
#' [beta_quadratic.DistribPenalty()] describes the points and the multivariate
#' case.
#'
#' # The meaning of `TRUE`
#'
#' `TRUE` does not imply that the marginal derivatives are available. A lasso
#' gives `TRUE`, because its log-density is linear in the response away from
#' the kink, and [penalty_dhessian()] still signals an error for it: the third
#' derivative is zero where it exists and undefined at the kink, which is where
#' a selecting fit places its coefficients.
#'
#' @param pen A [penalty()] object.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`. The
#'   generic checks it as [penalty_dhessian()] does, so a missing
#'   hyperparameter or a value outside its open bounds is an error on every
#'   branch. The separable branch evaluates the parent at it, and the other
#'   branches return a constant without reading its values.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A single logical. The base class returns `FALSE`, so a penalty
#'   class without a method of its own is treated as needing the third
#'   derivative.
#'
#' @seealso [is_quadratic()] for the property of the construction,
#'   [penalty_dhessian_beta()] for the derivative that is zero when this is
#'   `TRUE`, [penalty_hessian()] for the Hessian that is then constant in the
#'   coefficients.
#'
#' @examples
#' # Two of the three branches that give TRUE for both predicates.
#' beta_quadratic(quadratic_penalty(diag(3)), list(lambda = 1))
#' beta_quadratic(structured_penalty(
#'   parameters7::log_cholesky(2)),
#'   list(log_L1 = 0, log_L2 = 0, L2.1 = 0))
#'
#' # A Gaussian prior is quadratic in beta and is not a quadratic penalty.
#' g <- distrib_penalty(
#'   distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
#'   n_coef = 3)
#' c(beta_quadratic(g, list(sigma = 1)), is_quadratic(g))
#'
#' # A Student t prior is neither.
#' beta_quadratic(heavy_penalty(n_coef = 3), list(sigma = 1, nu = 4))
#'
#' # And TRUE does not mean the marginal derivatives exist.
#' beta_quadratic(lasso_penalty(n_coef = 3), list(lambda = 1))
#' try(penalty_dhessian(lasso_penalty(n_coef = 3), c(1, 0, -1),
#'                      list(lambda = 1)))
#'
#' @export
beta_quadratic <- S7::new_generic("beta_quadratic", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })


# --- the base class rejects -------------------------------------------------

#' @title Base Class Methods of the Marginal Generics
#' @name penalty_dhessian.penalty
#'
#' @description
#' The methods of the base class for the four generics that a marginal
#' criterion requires beyond the second order. `penalty_dhessian()`,
#' `penalty_d2hessian()` and `penalty_dcross()` each signal an error that
#' names the penalty and the generic; `beta_quadratic()` returns `FALSE`.
#' [scad_penalty()] and [mcp_penalty()] register none of the four and reach
#' these methods.
#'
#' @details
#' A marginal criterion based on a Laplace approximation requires derivatives
#' of the penalty beyond the second order. For a penalty that does not supply
#' them, these methods signal an error, so the hyperparameters of such a
#' penalty cannot be estimated by a marginal criterion. The message of
#' `penalty_dhessian()` states that consequence; the messages of
#' `penalty_d2hessian()` and `penalty_dcross()` name only the penalty and the
#' generic.
#'
#' The `FALSE` of `beta_quadratic()` is the conservative default: a penalty
#' class without a method of its own is treated as having a third
#' \eqn{\beta}-derivative.
#'
#' @param pen A [penalty()] object of a class that registers none of these
#'   methods.
#' @param beta A numeric vector of coefficients. Unused.
#' @param theta A named list of hyperparameter values. Unused.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_dhessian()`, `penalty_d2hessian()` and `penalty_dcross()`
#'   signal an error. `beta_quadratic()` returns the single logical `FALSE`.
#'
#' @seealso [penalty_dhessian()] and its siblings for the generics,
#'   [penalty_dhessian.QuadraticPenalty()] for a branch that supplies them.
#' @keywords internal
S7::method(penalty_dhessian, penalty) <- function(pen, beta, theta, ...) {
  stop(sprintf(paste0("'%s' does not supply penalty_dhessian(), so its\n",
                      "  hyperparameters cannot be estimated by a marginal",
                      " criterion."),
               pen@penalty_name), call. = FALSE)
}

#' @rdname penalty_dhessian.penalty
#' @name penalty_d2hessian.penalty
#' @keywords internal
S7::method(penalty_d2hessian, penalty) <- function(pen, beta, theta, ...) {
  stop(sprintf("'%s' does not supply penalty_d2hessian().",
               pen@penalty_name), call. = FALSE)
}

#' @rdname penalty_dhessian.penalty
#' @name penalty_dcross.penalty
#' @keywords internal
S7::method(penalty_dcross, penalty) <- function(pen, beta, theta, ...) {
  stop(sprintf("'%s' does not supply penalty_dcross().",
               pen@penalty_name), call. = FALSE)
}

#' @rdname penalty_dhessian.penalty
#' @name beta_quadratic.penalty
#' @keywords internal
S7::method(beta_quadratic, penalty) <- function(pen, theta, ...) FALSE


# --- quadratic --------------------------------------------------------------

#' @title Marginal Derivatives of a Quadratic Penalty
#' @name penalty_dhessian.QuadraticPenalty
#'
#' @description
#' The methods of the quadratic branch for the four generics that a marginal
#' criterion requires beyond the second order. \eqn{S = \lambda D'PD} is
#' linear in \eqn{\lambda} and free of the coefficients, so
#' `penalty_dhessian()` returns \eqn{D'PD}, `penalty_d2hessian()` and
#' `penalty_dcross()` return zeros, and `beta_quadratic()` returns `TRUE`.
#'
#' @details
#' `penalty_dhessian()` returns the stored \eqn{D'PD} **in the storage that
#' the penalty keeps**. A base matrix is passed through `unclass()`, which
#' removes the attributes that it may carry, and an S4 matrix is returned
#' unchanged. A penalty built with `blocks > 1` therefore gives a
#' `dgCMatrix`, and a consumer that needs a base matrix coerces it.
#'
#' The zeros are exact. The Hessian is \eqn{\lambda} times a constant matrix,
#' so its second derivative in \eqn{\lambda} vanishes identically, and so does
#' the third derivative of the value in \eqn{\beta} and two hyperparameters.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Unused, because none
#'   of these quantities depends on the coefficients.
#' @param theta A named list holding `lambda`. Unused, because none of these
#'   quantities depends on \eqn{\lambda}.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_dhessian()` a list of one matrix named `lambda`.
#'   `penalty_d2hessian()` a list of one zero matrix named `lambda_lambda`.
#'   `penalty_dcross()` a list of one zero vector named `lambda_lambda`.
#'   `beta_quadratic()` the single logical `TRUE`.
#'
#' @examples
#' pen <- quadratic_penalty(crossprod(diff(diag(3))))
#' b <- c(1, -0.5, 0.3)
#'
#' # The derivative is D'PD, so it does not change with lambda.
#' max(abs(penalty_dhessian(pen, b, list(lambda = 2))$lambda -
#'           penalty_dhessian(pen, b, list(lambda = 99))$lambda))
#'
#' # And everything above first order is exactly zero.
#' max(abs(unlist(penalty_d2hessian(pen, b, list(lambda = 2)))))
#' max(abs(unlist(penalty_dcross(pen, b, list(lambda = 2)))))
#' beta_quadratic(pen, list(lambda = 2))
#'
#' @seealso [penalty_dhessian()] and its siblings for the generics,
#'   [quadratic_penalty()] for the branch, [zero_pairs()] for the zeros,
#'   [penalty_dhessian_beta.QuadraticPenalty()] and
#'   [penalty_d2hessian_beta.QuadraticPenalty()] for the derivatives of the
#'   Hessian in the coefficients.
#' @keywords internal
S7::method(penalty_dhessian, QuadraticPenalty) <- function(pen, beta, theta,
                                                           ...) {
  # `unclass()` strips the attributes a base matrix may carry and does
  # nothing to an S4 one, which is what a blocked penalty stores: the
  # derivative is the matrix itself, in whatever storage the penalty keeps
  # it, and a caller that needs a base matrix coerces where the two meet.
  list(lambda = if (isS4(pen@DPD)) pen@DPD else unclass(pen@DPD))
}

#' @rdname penalty_dhessian.QuadraticPenalty
#' @name penalty_d2hessian.QuadraticPenalty
#' @keywords internal
S7::method(penalty_d2hessian, QuadraticPenalty) <- function(pen, beta, theta,
                                                            ...) {
  zero_pairs(pen, matrix(0, as.integer(pen@n_coef), as.integer(pen@n_coef)))
}

#' @rdname penalty_dhessian.QuadraticPenalty
#' @name penalty_dcross.QuadraticPenalty
#' @keywords internal
S7::method(penalty_dcross, QuadraticPenalty) <- function(pen, beta, theta,
                                                         ...) {
  zero_pairs(pen, numeric(as.integer(pen@n_coef)))
}

#' @rdname penalty_dhessian.QuadraticPenalty
#' @name beta_quadratic.QuadraticPenalty
#' @keywords internal
S7::method(beta_quadratic, QuadraticPenalty) <- function(pen, theta, ...) TRUE


# --- additive ---------------------------------------------------------------

#' @title Marginal Derivatives of an Additive Penalty
#' @name penalty_dhessian.AdditivePenalty
#'
#' @description
#' The methods of the additive branch for the four generics that a marginal
#' criterion requires beyond the second order. With \eqn{S_k = D'P_kD}, the
#' Hessian \eqn{S = \sum_k \lambda_k S_k} is linear in the smoothing
#' parameters and free of the coefficients, so `penalty_dhessian()` returns
#' the components \eqn{S_k} themselves, `penalty_d2hessian()` and
#' `penalty_dcross()` return zeros, and `beta_quadratic()` returns `TRUE`.
#'
#' @details
#' Every entry is exact, and none is computed at call time. The components
#' were fixed at construction and are stored already carried through the map,
#' so the first derivative is a lookup of `pen@mats` and the higher ones are
#' zeros of the right shape. Under the identity map \eqn{S_k = P_k}.
#'
#' [is_quadratic()] returns `TRUE` for this branch, as it does for the
#' quadratic and structured branches.
#'
#' @param pen An [AdditivePenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Unused.
#' @param theta A named list of `lambda1`, `lambda2`, ... Unused.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_dhessian()` a list of one base matrix per component, named
#'   by `pen@params`.
#'   `penalty_d2hessian()` a list of zero matrices, one per unordered pair.
#'   `penalty_dcross()` a list of zero vectors, one per unordered pair.
#'   `beta_quadratic()` the single logical `TRUE`.
#'
#' @examples
#' add <- additive_penalty(list(crossprod(diff(diag(4))), diag(4)))
#' b <- c(1, -0.5, 0.3, 0.2)
#' th <- list(lambda1 = 2, lambda2 = 0.5)
#'
#' # The derivative in each parameter is that parameter's own component.
#' d <- penalty_dhessian(add, b, th)
#' names(d)
#' max(abs(d$lambda2 - diag(4)))
#'
#' # Everything above first order is zero, over all three pairs.
#' names(penalty_d2hessian(add, b, th))
#' max(abs(unlist(penalty_d2hessian(add, b, th))))
#'
#' @seealso [penalty_dhessian()] and its siblings for the generics,
#'   [additive_penalty()] for the branch,
#'   [penalty_dhessian_beta.AdditivePenalty()] and
#'   [penalty_d2hessian_beta.AdditivePenalty()] for the derivatives of the
#'   Hessian in the coefficients.
#' @keywords internal
S7::method(penalty_dhessian, AdditivePenalty) <- function(pen, beta, theta,
                                                          ...) {
  stats::setNames(lapply(pen@mats, function(P)
    if (methods::is(P, "Matrix")) P else unclass(P)), pen@params)
}

#' @rdname penalty_dhessian.AdditivePenalty
#' @name penalty_d2hessian.AdditivePenalty
#' @keywords internal
S7::method(penalty_d2hessian, AdditivePenalty) <- function(pen, beta, theta,
                                                           ...) {
  zero_pairs(pen, additive_zero(pen))
}

#' @rdname penalty_dhessian.AdditivePenalty
#' @name penalty_dcross.AdditivePenalty
#' @keywords internal
S7::method(penalty_dcross, AdditivePenalty) <- function(pen, beta, theta,
                                                        ...) {
  zero_pairs(pen, numeric(as.integer(pen@n_coef)))
}

#' @rdname penalty_dhessian.AdditivePenalty
#' @name beta_quadratic.AdditivePenalty
#' @keywords internal
S7::method(beta_quadratic, AdditivePenalty) <- function(pen, theta, ...) TRUE


# --- structured -------------------------------------------------------------

#' @title Marginal Derivatives of a Structured Penalty
#' @name penalty_dhessian.StructuredPenalty
#'
#' @description
#' The methods of the structured branch for the four generics that a marginal
#' criterion requires beyond the second order. \eqn{S = \Omega(\theta)} is
#' the matrix parameter itself, so `penalty_dhessian()` and
#' `penalty_d2hessian()` are the derivatives `param_d1` and `param_d2` of the
#' structure, `penalty_dcross()` is the second of those applied to the
#' coefficients, and `beta_quadratic()` is `TRUE`.
#'
#' @details
#' Nothing is differentiated here. \pkg{parameters7} supplies both derivative
#' arrays in closed form for the structures that it ships, and through its
#' numerical fallback for a structure that does not implement them; this
#' branch unwraps and re-keys the arrays. The free values in `theta` are
#' passed to the structure as they are, in `pen@params` order.
#'
#' Since the value is quadratic in the coefficients, the mixed third derivative
#' is \eqn{(\partial^2\Omega/\partial\theta_m\partial\theta_l)\beta}, so
#' `penalty_dcross()` calls `penalty_d2hessian()` and multiplies.
#'
#' # The keying
#'
#' `param_d2` is a **flat** list whose keys join the free names with a
#' separator of its own, so an entry is located through
#' [parameters7::param_tuple_indices()], which enumerates the tuples in the
#' same order as the components. The method does not parse the keys, because
#' a free name may itself contain the separator. If a pair has no component,
#' the method signals an error.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Read by
#'   `penalty_dcross()`; the other three methods ignore it.
#' @param theta A named list of the free values of the structure.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_dhessian()` a list of one matrix per free value, named by
#'   `pen@params` and with dimnames stripped.
#'   `penalty_d2hessian()` a list of one matrix per unordered pair, keyed
#'   diagonals first.
#'   `penalty_dcross()` a list of one numeric vector of length `pen@n_coef` per
#'   unordered pair, keyed the same way.
#'   `beta_quadratic()` the single logical `TRUE`.
#'
#' @examples
#' s <- structured_penalty(parameters7::log_cholesky(2))
#' th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
#' b <- c(1, -0.5)
#'
#' # Three free values, so three first derivatives and six pairs.
#' names(penalty_dhessian(s, b, th))
#' names(penalty_d2hessian(s, b, th))
#'
#' # Unlike the quadratic branch, the second order does not vanish.
#' penalty_d2hessian(s, b, th)$log_L1_log_L1
#'
#' # And the mixed block is that matrix applied to the coefficients.
#' all.equal(penalty_dcross(s, b, th)$log_L1_log_L1,
#'           drop(penalty_d2hessian(s, b, th)$log_L1_log_L1 %*% b))
#'
#' @seealso [penalty_dhessian()] and its siblings for the generics,
#'   [structured_penalty()] for the branch,
#'   [parameters7::param_d1()] for the arrays that these methods read,
#'   [penalty_dhessian_beta.StructuredPenalty()] and
#'   [penalty_d2hessian_beta.StructuredPenalty()] for the derivatives of the
#'   Hessian in the coefficients.
#' @keywords internal
S7::method(penalty_dhessian, StructuredPenalty) <- function(pen, beta, theta,
                                                            ...) {
  d1 <- parameters7::param_d1(pen@structure, struct_eta(pen, theta))
  stats::setNames(lapply(seq_along(pen@params), function(m)
    unname(unclass(d1[[m]]))), pen@params)
}

#' @rdname penalty_dhessian.StructuredPenalty
#' @name penalty_d2hessian.StructuredPenalty
#' @keywords internal
S7::method(penalty_d2hessian, StructuredPenalty) <- function(pen, beta, theta,
                                                             ...) {
  d2 <- parameters7::param_d2(pen@structure, struct_eta(pen, theta))
  # param_d2 is a FLAT list, not a nested one, and its keys join the free
  # names with a separator of its own. The entry is located through
  # param_tuple_indices(), which enumerates the tuples in exactly the order
  # the names are in, so nothing here depends on the spelling of a key.
  idx <- parameters7::param_tuple_indices(pen@structure, 2L)
  prs <- ptheta_pairs(pen@params)
  stats::setNames(lapply(names(prs), function(nm) {
    ij <- sort(match(prs[[nm]], pen@params))
    hit <- which(vapply(idx, function(t) identical(sort(t), ij), logical(1)))
    if (!length(hit)) {
      stop(sprintf("No second-derivative component for '%s' and '%s'.",
                   prs[[nm]][1L], prs[[nm]][2L]), call. = FALSE)
    }
    unname(unclass(d2[[hit[1L]]]))
  }), names(prs))
}

#' @rdname penalty_dhessian.StructuredPenalty
#' @name penalty_dcross.StructuredPenalty
#' @keywords internal
S7::method(penalty_dcross, StructuredPenalty) <- function(pen, beta, theta,
                                                          ...) {
  d2 <- penalty_d2hessian(pen, beta, theta)
  lapply(d2, function(M) as.numeric(M %*% beta))
}

#' @rdname penalty_dhessian.StructuredPenalty
#' @name beta_quadratic.StructuredPenalty
#' @keywords internal
S7::method(beta_quadratic, StructuredPenalty) <- function(pen, theta, ...) TRUE


# --- separable --------------------------------------------------------------

#' @title Marginal Derivatives of a Separable Penalty
#' @name penalty_dhessian.DistribPenalty
#'
#' @description
#' The methods of the separable branch for the four generics that a marginal
#' criterion requires beyond the second order. With
#' \eqn{\rho = -\sum_j \log f((D\beta)_j;\theta)} the Hessian is
#' \eqn{-D'\mathrm{diag}(\ell^{(yy)})D}, so its \eqn{\theta}-derivatives are
#' the higher components of the parent carried through the same map.
#'
#' @details
#' # The component of the parent behind each method
#'
#' | generic | component of the parent | result |
#' |---|---|---|
#' | `penalty_dhessian()` | [distributions7::distrib_cross2_y()] | \eqn{-D'\mathrm{diag}(\partial^3\ell/\partial y^2\partial\theta_m)D} |
#' | `penalty_d2hessian()` | [distributions7::distrib_hess_y_hess()] | \eqn{-D'\mathrm{diag}(\partial^4\ell/\partial y^2\partial\theta_m\partial\theta_l)D} |
#' | `penalty_dcross()` | [distributions7::distrib_grad_y_hess()] | \eqn{-D'(\partial^3\ell/\partial y\,\partial\theta_m\partial\theta_l)} |
#'
#' Each component is evaluated at \eqn{D\beta} and at `theta` on the parameter
#' scale of the parent, so the derivatives are in \eqn{\theta} itself and no
#' link is applied.
#'
#' Nothing is differentiated here. A parent with closed forms for these
#' components, such as the Gaussian and therefore every Gaussian random
#' effect, makes this branch exact. A parent without them takes the fallback
#' of \pkg{distributions7}, which applies one stencil to the highest-order
#' quantity that the family implements itself, usually its analytic response
#' Hessian.
#'
#' A multivariate parent is read blockwise: each component is assembled into a
#' block-diagonal matrix by [dp_blockdiag()] and carried through the map by
#' [map_quad_full()], and a univariate one goes through [map_quad()].
#'
#' # Kinked parents
#'
#' A penalty whose parent declares a kink, which is the lasso and the elastic
#' net, signals an error in `penalty_dhessian()`, `penalty_d2hessian()` and
#' `penalty_dcross()`, with the penalty named (see [reject_kinked()]). The
#' third derivative does not exist at the kink, and the mode at which a
#' marginal criterion is expanded is where a kinked penalty places
#' coefficients. `beta_quadratic()` signals no error and returns `TRUE` for
#' the lasso, whose log-density is linear in the response away from the kink.
#'
#' # The test of `beta_quadratic()`
#'
#' The log-density is quadratic in the response exactly when its second
#' derivative in the response does not depend on the response. The generic
#' has already checked `theta` (an invalid value is an error), and the method
#' evaluates [distributions7::distrib_hess_y()] of the parent at four points of
#' `pen@block` coordinates each, and tests with [all.equal()] at a relative
#' tolerance of `1e-12` whether the values agree across the points. A
#' multivariate parent returns one \eqn{p \times p} matrix per point, or a
#' single matrix when that matrix does not depend on the response, and the
#' comparison is made between the matrices of different points, not between
#' the entries of one matrix. A single value or a single matrix, which carries
#' no dimension for the points, gives `TRUE`. If the parent signals an error,
#' or returns a number of values that is not a multiple of `pen@block^2`, the
#' method returns `FALSE`.
#'
#' @param pen A [DistribPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Read by all three
#'   derivative methods, because the components of the parent depend on it.
#' @param theta A named list of the free parameters of the parent.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_dhessian()` a list of one base matrix of side `pen@n_coef`
#'   per hyperparameter.
#'   `penalty_d2hessian()` a list of one such matrix per unordered pair, keyed
#'   diagonals first.
#'   `penalty_dcross()` a list of one numeric vector of length `pen@n_coef` per
#'   unordered pair.
#'   `beta_quadratic()` a single logical.
#'
#' @examples
#' b <- c(1, -0.5, 0.3)
#'
#' # A Gaussian prior at zero: the Hessian is 1/sigma^2 on the diagonal, so
#' # its derivative in sigma is -2/sigma^3.
#' g <- distrib_penalty(
#'   distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
#'   n_coef = 3)
#' diag(penalty_dhessian(g, b, list(sigma = 1.2))$sigma)
#' rep(-2 / 1.2^3, 3)
#'
#' # A Student t prior is not quadratic in the coefficients, so its
#' # derivative changes with them.
#' h <- heavy_penalty(n_coef = 3)
#' beta_quadratic(h, list(sigma = 1, nu = 4))
#' diag(penalty_dhessian(h, b, list(sigma = 1, nu = 4))$nu)
#'
#' # A kinked parent has no third derivative, and the call is an error.
#' try(penalty_dhessian(lasso_penalty(n_coef = 3), b, list(lambda = 1)))
#'
#' @seealso [penalty_dhessian()] and its siblings for the generics,
#'   [distrib_penalty()] for the branch, [reject_kinked()] for the error on a
#'   kinked parent, [carry_pairs()] for the re-keying,
#'   [penalty_dhessian_beta.DistribPenalty()] and
#'   [penalty_d2hessian_beta.DistribPenalty()] for the derivatives of the
#'   Hessian in the coefficients.
#' @keywords internal
S7::method(penalty_dhessian, DistribPenalty) <- function(pen, beta, theta,
                                                         ...) {
  reject_kinked(pen, "penalty_dhessian")
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    c2 <- distributions7::distrib_cross2_y(pen@parent, a, theta)
    return(stats::setNames(lapply(pen@params, function(m)
      -map_quad_full(pen, dp_blockdiag(pen, c2[[m]], nrow(a)))), pen@params))
  }
  c2 <- distributions7::distrib_cross2_y(pen@parent, t, theta)
  stats::setNames(lapply(pen@params, function(m)
    -map_quad(pen, c2[[m]] + 0 * t)), pen@params)
}

#' @rdname penalty_dhessian.DistribPenalty
#' @name penalty_d2hessian.DistribPenalty
#' @keywords internal
S7::method(penalty_d2hessian, DistribPenalty) <- function(pen, beta, theta,
                                                          ...) {
  reject_kinked(pen, "penalty_d2hessian")
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    h <- distributions7::distrib_hess_y_hess(pen@parent, a, theta)
    return(carry_pairs(pen, h, function(v)
      -map_quad_full(pen, dp_blockdiag(pen, v, nrow(a)))))
  }
  h <- distributions7::distrib_hess_y_hess(pen@parent, t, theta)
  carry_pairs(pen, h, function(v) -map_quad(pen, v + 0 * t))
}

#' @rdname penalty_dhessian.DistribPenalty
#' @name penalty_dcross.DistribPenalty
#' @keywords internal
S7::method(penalty_dcross, DistribPenalty) <- function(pen, beta, theta, ...) {
  reject_kinked(pen, "penalty_dcross")
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    g <- distributions7::distrib_grad_y_hess(pen@parent, a, theta)
    return(carry_pairs(pen, g, function(v)
      -map_back(pen, dp_flat(pen, v + 0 * a))))
  }
  g <- distributions7::distrib_grad_y_hess(pen@parent, t, theta)
  carry_pairs(pen, g, function(v) -map_back(pen, v + 0 * t))
}

#' @rdname penalty_dhessian.DistribPenalty
#' @name beta_quadratic.DistribPenalty
#' @keywords internal
S7::method(beta_quadratic, DistribPenalty) <- function(pen, theta, ...) {
  # the log-density is quadratic in the response exactly when its second
  # derivative there does not depend on it. The question is put that way, and
  # not to a third derivative, because the second is analytic for almost every
  # family while the third is often a difference, whose noise no threshold
  # separates from a true zero.
  #
  # ⚠️ THE PROBE IS SHAPED BY dp_arg(), and the comparison is between
  # OBSERVATIONS rather than between the entries of one reading. A
  # multivariate parent answers with one p by p matrix per observation, or
  # with a single one where that matrix does not move -- which is itself the
  # proof of quadraticity -- so asking whether every entry equals the first
  # compares an off-diagonal with a diagonal and is false for any prior over
  # more than one coordinate. Measured on the multivariate gaussian a
  # covariance class declares: hess_y is the constant -Sigma^-1, whose
  # entries read -1, 0, 0, -1, and the old form returned FALSE where the
  # prior is exactly quadratic. A slice is k = block^2 numbers, one for a
  # univariate parent, so the two cases are one expression.
  # FOUR PROBE POINTS, of pen@block coordinates each. A p-variate parent
  # reads a point as p numbers, so four VALUES laid out with p columns give
  # two rows at p = 2 or 3 and ONE from p = 4 up -- and a single slice takes
  # the branch below, which answers TRUE without comparing anything.
  # Measured on a multivariate Student t prior, whose hess_y does move with
  # the response: four values answer FALSE at p = 2 and 3 and TRUE at 4 and
  # 5, where the truth is FALSE at every p. The offset is zero at p = 1, so
  # a univariate parent is probed at exactly the four values it always was.
  # theta is validated by the generic, outside the tryCatch, which is there
  # for a parent that cannot evaluate its response curvature and not for an
  # invalid argument
  p <- as.integer(pen@block)
  t <- dp_arg(pen, rep(c(-1.73, -0.29, 0.61, 2.04), each = p) +
                     0.37 * (rep(seq_len(p), 4L) - 1L))
  h <- tryCatch(distributions7::distrib_hess_y(pen@parent, t, theta),
                error = function(e) NULL)
  if (is.null(h) || !length(h)) return(FALSE)
  k <- as.integer(pen@block)^2L
  if (length(h) %% k) return(FALSE)
  n <- length(h) %/% k
  # a reading carrying no observation dimension does not move with the
  # response, which is the answer rather than a case to compare
  if (n <= 1L) return(TRUE)
  v <- as.numeric(h)
  isTRUE(all.equal(v, rep(v[seq_len(k)], n), tolerance = 1e-12))
}


#' Carry a Parent's Paired Components Into Coefficient Space
#'
#' @description
#' Re-keys a \pkg{distributions7} component keyed by parameter pair into this
#' package's own keys, and places each through the penalty's map.
#'
#' @details
#' The two enumerations of pairs are built the same way from the same names, so
#' a key of one is a key of the other. Each pair is looked up by name, first in
#' `pen@params` order and then in the reversed order, and not by position or
#' by splitting a key, since a hyperparameter whose own name contains the
#' separator could not be split back into its two names. If neither name is
#' present, the function signals an error that names the pair.
#'
#' @param pen A [penalty()] object.
#' @param comp The components of the parent, keyed by parameter pair.
#' @param carry A function that places one component into coefficient space.
#'
#' @return A named list keyed by hyperparameter pair, in the order of
#'   [ptheta_pairs()].
#'
#' @keywords internal
carry_pairs <- function(pen, comp, carry) {
  prs <- ptheta_pairs(pen@params)
  keys <- names(comp)
  stats::setNames(lapply(names(prs), function(nm) {
    pr <- prs[[nm]]
    key <- paste(pr[1], pr[2], sep = "_")
    if (!key %in% keys) key <- paste(pr[2], pr[1], sep = "_")
    if (!key %in% keys) {
      stop(sprintf("The parent has no component for '%s' and '%s'.",
                   pr[1], pr[2]), call. = FALSE)
    }
    carry(comp[[key]])
  }), names(prs))
}


#' A Zero Entry for Every Hyperparameter Pair
#'
#' @description
#' Returns a list holding the same zero object under every hyperparameter-pair
#' key, which is the value of [penalty_d2hessian()] and [penalty_dcross()] for
#' a penalty whose Hessian is linear in its hyperparameters and free of the
#' coefficients.
#'
#' @details
#' The quadratic and additive branches use it for both
#' [penalty_d2hessian()] and [penalty_dcross()], the two differing only in the
#' shape of the zero.
#'
#' @param pen A [penalty()] object.
#' @param z The zero object: a `pen@n_coef` by `pen@n_coef` matrix for
#'   [penalty_d2hessian()], a numeric vector of that length for
#'   [penalty_dcross()].
#'
#' @return A named list of \eqn{p(p+1)/2} entries, keyed as
#'   [ptheta_pairs()] keys them, each holding `z`.
#'
#' @seealso [penalty_d2hessian()], [penalty_dcross()], [ptheta_pairs()]
#'
#' @keywords internal
zero_pairs <- function(pen, z) {
  nm <- names(ptheta_pairs(pen@params))
  stats::setNames(rep(list(z), length(nm)), nm)
}


#' Reject a Kinked Penalty
#'
#' @description
#' Signals an error naming the penalty and the generic when `pen@kinks` holds
#' at least one finite value, and returns invisibly otherwise. It is called at
#' the head of the six methods of the separable branch that need a third or
#' higher derivative: `penalty_dhessian()`, `penalty_d2hessian()`,
#' `penalty_dcross()`, `penalty_dhessian_beta()`, `penalty_d2hessian_beta()`
#' and `penalty_dhessian_beta_theta()`.
#'
#' @details
#' A marginal criterion is a Laplace expansion at the penalized mode, and a
#' kinked penalty places coefficients exactly at the kink, where the third
#' derivative does not exist. The function therefore signals an error and
#' does not return a one-sided derivative.
#'
#' The kinks are read directly from `pen@kinks` and not through
#' [penalty_kinks()], because that generic takes a `theta` and this function
#' receives only the penalty and the name of the generic. A penalty without a
#' `kinks` property is treated as having no kink; every caller in the package
#' is a method of [DistribPenalty()], which has the property.
#'
#' @param pen A [penalty()] object.
#' @param what The name of the generic, a single string, inserted in the
#'   message followed by `()`.
#'
#' @return `NULL`, invisibly, when there is no kink.
#'
#' @section Errors:
#' `'<penalty>' has a kink, so <what>() does not exist there and its
#' hyperparameters cannot be estimated by a marginal criterion.`
#'
#' @seealso [penalty_dhessian.DistribPenalty()], [penalty_kinks()]
#'
#' @keywords internal
reject_kinked <- function(pen, what) {
  # the property is read directly rather than through penalty_kinks(), which
  # takes a theta this is asked before there is one to pass
  k <- if ("kinks" %in% S7::prop_names(pen)) pen@kinks else numeric(0)
  if (length(k) && any(is.finite(unlist(k)))) {
    stop(sprintf(paste0("'%s' has a kink, so %s() does not exist there and\n",
                        "  its hyperparameters cannot be estimated by a",
                        " marginal criterion."),
                 pen@penalty_name, what), call. = FALSE)
  }
  invisible(NULL)
}


# --- the derivative of the coefficient Hessian in the coefficients -----------

#' The Derivative of the Coefficient Hessian in the Coefficients
#'
#' @description
#' Returns \eqn{\partial S/\partial\beta\,[v] = \sum_c v_c\,
#' \partial^3\rho/\partial\beta\,\partial\beta'\,\partial\beta_c}, the
#' derivative of the coefficient Hessian along a direction \eqn{v}, as one
#' matrix of side `pen@n_coef`. It is zero for every penalty whose Hessian
#' does not depend on the coefficients, apart from a kinked separable penalty,
#' for which the generic signals an error.
#'
#' @details
#' # Role in a marginal criterion
#'
#' A marginal criterion reads \eqn{\log\lvert H + S\rvert} at the penalized
#' mode, and the mode moves with the hyperparameters. When \eqn{S} does not
#' depend on the coefficients, that movement reaches the determinant through
#' \eqn{H} alone. When it does, as for a heavy-tailed prior on a random effect,
#' the movement reaches the determinant through \eqn{S} as well, and the
#' gradient contains \eqn{\mathrm{tr}(M\,\partial S/\partial\beta[v])} with
#' \eqn{v} the derivative of the mode. A prediction-error criterion contains
#' the same matrix inside the trace that gives its degrees of freedom.
#'
#' The method contracts the direction and does not return the array of order
#' three, which has \eqn{q^3} entries, because the criteria read that array
#' only along a direction.
#'
#' # The value on each branch
#'
#' | branch | \eqn{\partial S/\partial\beta[v]} |
#' |---|---|
#' | [quadratic_penalty()], [additive_penalty()], [structured_penalty()] | zero |
#' | [distrib_penalty()] whose parent is quadratic in the argument | zero |
#' | [distrib_penalty()] with a univariate parent otherwise | \eqn{-D'\mathrm{diag}(\ell^{(yyy)}(D\beta)\odot Dv)\,D} |
#' | [distrib_penalty()] with a multivariate parent otherwise | \eqn{-D'\mathrm{blockdiag}(T_i)D}, \eqn{(T_i)_{ab} = \sum_c \ell^{(y_ay_by_c)}(b_i)(Dv)_{ic}} |
#' | a kinked parent, [scad_penalty()], [mcp_penalty()] | an error |
#'
#' The univariate row follows from \eqn{S = -D'\mathrm{diag}(\ell^{(yy)}(D\beta))D}:
#' differentiating \eqn{\ell^{(yy)}((D\beta)_j)} in \eqn{\beta_c} gives
#' \eqn{\ell^{(yyy)}((D\beta)_j) D_{jc}}, and summing against \eqn{v_c} gives
#' the factor \eqn{(Dv)_j}. The third response derivative is read from
#' [distributions7::distrib_deriv3_y()], which is closed form for every family
#' that \pkg{distributions7} ships, the Student t among them, and one stencil
#' on the log-density for a family defined by its density alone.
#'
#' In the multivariate row \eqn{b_i} is block \eqn{i} of \eqn{D\beta} and
#' \eqn{(Dv)_{ic}} is coordinate \eqn{c} of block \eqn{i} of \eqn{Dv}. A
#' multivariate parent that is not quadratic, such as a multivariate Student t
#' prior, is read through [distributions7::distrib_deriv3_y()], which returns
#' one \eqn{p \times p \times p} array per block; [dp_contract()] contracts its
#' last index against the coordinates of the block in \eqn{Dv}. If the parent
#' cannot supply the array, the method signals an error that names the
#' penalty (see [dp_parent_or_reject()]). The separable method calls
#' [beta_quadratic()] first, which evaluates the response Hessian of the
#' parent, and returns zero when it gives `TRUE`. A Gaussian prior of any
#' dimension therefore gives zero, and the third response derivative of the
#' parent is not computed.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`; the generic signals an
#'   error for any other length.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, checked by the generic as in
#'   [penalty_dhessian()].
#' @param v A numeric vector of length `pen@n_coef`, the direction; the
#'   generic signals an error for any other length.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A square base matrix of side `pen@n_coef`.
#'
#' @seealso [penalty_hessian()] for the quantity differentiated,
#'   [beta_quadratic()] for the predicate under which it is zero,
#'   [penalty_dhessian()] for the derivative in the hyperparameters,
#'   [penalty_d2hessian_beta()] for the derivatives of this quantity.
#'
#' @examples
#' b <- c(1, -0.5, 0.3)
#' v <- c(0.2, 0.1, -0.4)
#'
#' # Zero for a quadratic penalty, whose Hessian does not change with beta.
#' penalty_dhessian_beta(quadratic_penalty(diag(3)), b, list(lambda = 2), v)
#'
#' # The Hessian of a Student t prior changes with beta, and this is its
#' # derivative along v, compared with a central difference.
#' h <- heavy_penalty(n_coef = 3)
#' th <- list(sigma = 1, nu = 4)
#' D <- penalty_dhessian_beta(h, b, th, v)
#' eps <- 1e-6
#' num <- (penalty_hessian(h, b + eps * v, th) -
#'         penalty_hessian(h, b - eps * v, th)) / (2 * eps)
#' max(abs(D - num))
#'
#' # A kinked parent has no such derivative at the kink, and the call is an
#' # error.
#' try(penalty_dhessian_beta(lasso_penalty(n_coef = 3), b, list(lambda = 1), v))
#'
#' @export
penalty_dhessian_beta <- S7::new_generic("penalty_dhessian_beta", "pen",
  function(pen, beta, theta, v, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    v <- as.numeric(v)
    if (length(v) != length(beta)) {
      stop(sprintf("'v' has length %d where the coefficients have %d.",
                   length(v), length(beta)), call. = FALSE)
    }
    S7::S7_dispatch()
  })

#' @title Methods of penalty_dhessian_beta()
#' @name penalty_dhessian_beta.penalty
#'
#' @description
#' The method of the base class signals an error that names the penalty; it is
#' the method that [scad_penalty()] and [mcp_penalty()] reach. The quadratic,
#' additive and structured branches return a zero matrix, because their
#' Hessian is free of the coefficients. The separable branch returns zero when
#' [beta_quadratic()] is `TRUE` for it, and
#' \eqn{-D'\mathrm{diag}(\ell^{(yyy)}(D\beta)\odot Dv)D} otherwise, with a
#' block-diagonal middle matrix for a multivariate parent. It signals an error
#' for a kinked parent and for a multivariate parent that cannot supply its
#' third response derivative.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of hyperparameter values.
#' @param v A numeric vector of length `pen@n_coef`, the direction.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A square base matrix of side `pen@n_coef`, or an error.
#'
#' @seealso [penalty_dhessian_beta()] for the generic.
#' @keywords internal
S7::method(penalty_dhessian_beta, penalty) <- function(pen, beta, theta, v,
                                                       ...) {
  stop(sprintf(paste0("'%s' does not supply penalty_dhessian_beta(), so the\n",
                      "  movement of its Hessian with the coefficients is not",
                      " available."),
               pen@penalty_name), call. = FALSE)
}

#' @rdname penalty_dhessian_beta.penalty
#' @name penalty_dhessian_beta.QuadraticPenalty
#' @keywords internal
S7::method(penalty_dhessian_beta, QuadraticPenalty) <- function(pen, beta,
                                                                theta, v, ...) {
  matrix(0, as.integer(pen@n_coef), as.integer(pen@n_coef))
}

#' @rdname penalty_dhessian_beta.penalty
#' @name penalty_dhessian_beta.AdditivePenalty
#' @keywords internal
S7::method(penalty_dhessian_beta, AdditivePenalty) <- function(pen, beta,
                                                               theta, v, ...) {
  additive_zero(pen)
}

#' @rdname penalty_dhessian_beta.penalty
#' @name penalty_dhessian_beta.StructuredPenalty
#' @keywords internal
S7::method(penalty_dhessian_beta, StructuredPenalty) <- function(pen, beta,
                                                                 theta, v,
                                                                 ...) {
  matrix(0, as.integer(pen@n_coef), as.integer(pen@n_coef))
}

#' @rdname penalty_dhessian_beta.penalty
#' @name penalty_dhessian_beta.DistribPenalty
#' @keywords internal
S7::method(penalty_dhessian_beta, DistribPenalty) <- function(pen, beta, theta,
                                                              v, ...) {
  reject_kinked(pen, "penalty_dhessian_beta")
  k <- as.integer(pen@n_coef)
  # a parent quadratic in its argument has a constant response curvature, so
  # the derivative is exactly zero and the parent is not asked for a third
  # derivative it may only have as a stencil, or not at all
  if (isTRUE(beta_quadratic(pen, theta))) return(matrix(0, k, k))
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    d3 <- dp_parent_or_reject(pen, "third", distributions7::distrib_deriv3_y,
                              a, theta)
    h <- dp_contract(d3, list(dp_arg(pen, map_apply(pen, v))))
    return(-map_quad_full(pen, dp_blockdiag(pen, h, nrow(a))))
  }
  d3 <- distributions7::distrib_deriv3_y(pen@parent, t, theta) + 0 * t
  -map_quad(pen, d3 * map_apply(pen, v))
}


#' Second Derivatives of the Coefficient Hessian in the Coefficients and the
#' Hyperparameters
#'
#' @description
#' The two derivatives of [penalty_dhessian_beta()] that the second derivative
#' of a marginal criterion contains. `penalty_d2hessian_beta()` is
#' \eqn{\partial^2 S/\partial\beta^2\,[v, w]}, the second derivative of the
#' Hessian in the coefficients contracted along two directions.
#' `penalty_dhessian_beta_theta()` is
#' \eqn{\partial^2 S/\partial\beta\,\partial\theta_m\,[v]}, the derivative of
#' \eqn{\partial S/\partial\beta\,[v]} in each hyperparameter.
#'
#' @details
#' Write \eqn{S_m = \partial S/\partial\theta_m}, \eqn{T[v] = \partial
#' S/\partial\beta\,[v]}, and \eqn{b_m} and \eqn{b_{ml}} for the first and
#' second derivatives of the mode in the hyperparameters. The matrix
#' \eqn{K_m = S_m + T[b_m]} enters the gradient of the criterion, and its
#' derivative in \eqn{\theta_l} has three parts that are present only for a
#' penalty whose Hessian depends on the coefficients: \eqn{T[b_{ml}]}, from
#' [penalty_dhessian_beta()]; the second derivative in the coefficients along
#' \eqn{b_m} and \eqn{b_l}, from `penalty_d2hessian_beta()`; and the two
#' matrices \eqn{\partial S_m/\partial\beta\,[b_l]} and
#' \eqn{\partial S_l/\partial\beta\,[b_m]}. Both of these come from
#' `penalty_dhessian_beta_theta()`, the first as its component \eqn{m} along
#' \eqn{b_l} and the second as its component \eqn{l} along \eqn{b_m}, because
#' by the symmetry of mixed partials in \eqn{\beta} and \eqn{\theta} the
#' derivative of \eqn{T[b_m]} in \eqn{\theta_l} equals
#' \eqn{\partial S_l/\partial\beta\,[b_m]}. The two matrices are in general
#' different.
#'
#' # The value on each branch
#'
#' | branch | \eqn{\partial^2 S/\partial\beta^2[v,w]} | \eqn{\partial^2 S/\partial\beta\partial\theta_m[v]} |
#' |---|---|---|
#' | quadratic, additive, structured | zero | zero for every hyperparameter |
#' | separable, parent quadratic in its argument | zero | zero for every hyperparameter |
#' | separable, univariate parent otherwise | \eqn{-D'\mathrm{diag}(\ell^{(yyyy)} \odot Dv \odot Dw)D} | \eqn{-D'\mathrm{diag}(\partial_{\theta_m}\ell^{(yyy)} \odot Dv)D} |
#' | separable, multivariate parent otherwise | the same per block, contracted by [dp_contract()] | the same per block |
#' | a kinked parent, [scad_penalty()], [mcp_penalty()] | an error | an error |
#'
#' The univariate rows follow from
#' \eqn{\partial S/\partial\beta[v] = -D'\mathrm{diag}(\ell^{(yyy)}(D\beta) \odot
#' Dv)D} by differentiating \eqn{\ell^{(yyy)}((D\beta)_j)} once more, in
#' \eqn{\beta} along \eqn{w} or in \eqn{\theta_m}. The fourth response
#' derivative is [distributions7::distrib_deriv4_y()], closed form for every
#' family that \pkg{distributions7} ships and one stencil for a family defined
#' by its density alone. The mixed one is
#' [distributions7::distrib_cross3_y()], closed form for the location families
#' listed on its page, the Student t among them, and one central difference of
#' the third response derivative for the other families. A multivariate
#' parent supplies the two as arrays per block, of \eqn{p^4} and \eqn{p^3}
#' entries, which [dp_contract()] contracts against the coordinates of the
#' block in the directions; the multivariate Student t supplies both. The
#' derivatives are taken in \eqn{\theta} on its own scale.
#'
#' The methods of the base class and of each branch are described on
#' [penalty_d2hessian_beta.penalty()].
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`; the generics signal an
#'   error for any other length.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, checked by the generics as in
#'   [penalty_dhessian()].
#' @param v,w Numeric vectors of length `pen@n_coef`, the directions; the
#'   generics signal an error for any other length.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return `penalty_d2hessian_beta()` a square base matrix of side
#'   `pen@n_coef`; `penalty_dhessian_beta_theta()` a list of such matrices, one
#'   per hyperparameter, keyed by `pen@params`.
#'
#' @seealso [penalty_dhessian_beta()] for the quantity differentiated,
#'   [penalty_dhessian()] for the derivative of the Hessian in the
#'   hyperparameters, [beta_quadratic()] for the predicate under which both
#'   are zero.
#'
#' @examples
#' h <- heavy_penalty(n_coef = 3)
#' b <- c(1, -0.5, 0.3)
#' v <- c(0.2, 0.1, -0.4)
#' w <- c(-0.3, 0.5, 0.1)
#' th <- list(sigma = 1, nu = 4)
#'
#' # The second derivative along two directions, against a difference of the
#' # first along one of them.
#' eps <- 1e-6
#' D2 <- penalty_d2hessian_beta(h, b, th, v, w)
#' num <- (penalty_dhessian_beta(h, b + eps * w, th, v) -
#'         penalty_dhessian_beta(h, b - eps * w, th, v)) / (2 * eps)
#' max(abs(D2 - num))
#'
#' # The movement in the degrees of freedom.
#' Dnu <- penalty_dhessian_beta_theta(h, b, th, v)$nu
#' num <- (penalty_dhessian_beta(h, b, list(sigma = 1, nu = 4 + eps), v) -
#'         penalty_dhessian_beta(h, b, list(sigma = 1, nu = 4 - eps), v)) /
#'   (2 * eps)
#' max(abs(Dnu - num))
#'
#' # Zero for a quadratic penalty.
#' penalty_d2hessian_beta(quadratic_penalty(diag(3)), b, list(lambda = 2), v, w)
#'
#' @export
penalty_d2hessian_beta <- S7::new_generic("penalty_d2hessian_beta", "pen",
  function(pen, beta, theta, v, w, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    v <- as.numeric(v)
    w <- as.numeric(w)
    if (length(v) != length(beta) || length(w) != length(beta)) {
      stop(sprintf(paste0("'v' and 'w' have lengths %d and %d where the",
                          " coefficients have %d."),
                   length(v), length(w), length(beta)), call. = FALSE)
    }
    S7::S7_dispatch()
  })

#' @rdname penalty_d2hessian_beta
#' @export
penalty_dhessian_beta_theta <- S7::new_generic("penalty_dhessian_beta_theta",
                                               "pen",
  function(pen, beta, theta, v, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    v <- as.numeric(v)
    if (length(v) != length(beta)) {
      stop(sprintf("'v' has length %d where the coefficients have %d.",
                   length(v), length(beta)), call. = FALSE)
    }
    S7::S7_dispatch()
  })

#' @title Methods of penalty_d2hessian_beta() and penalty_dhessian_beta_theta()
#' @name penalty_d2hessian_beta.penalty
#' @aliases penalty_d2hessian_beta.QuadraticPenalty
#'   penalty_d2hessian_beta.AdditivePenalty
#'   penalty_d2hessian_beta.StructuredPenalty
#'   penalty_dhessian_beta_theta.QuadraticPenalty
#'   penalty_dhessian_beta_theta.AdditivePenalty
#'   penalty_dhessian_beta_theta.StructuredPenalty
#'
#' @description
#' The methods of the base class for both generics signal an error that names
#' the penalty; they are the methods that [scad_penalty()] and [mcp_penalty()]
#' reach. The quadratic, additive and structured branches
#' ([QuadraticPenalty()], [AdditivePenalty()] and [StructuredPenalty()]) share
#' one method per generic, which returns a zero matrix of side `pen@n_coef`
#' from `penalty_d2hessian_beta()` and a list of such zero matrices, one per
#' hyperparameter, from `penalty_dhessian_beta_theta()`, because their Hessian
#' is free of the coefficients. The separable branch returns the same zeros
#' when [beta_quadratic()] is `TRUE` for it, and otherwise the closed forms
#' given on the page of [penalty_d2hessian_beta()], contracted per block for a
#' multivariate parent. It signals an error for a kinked parent and for a
#' multivariate parent that cannot supply its fourth or mixed third response
#' derivative.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A named list of hyperparameter values.
#' @param v,w Numeric vectors of length `pen@n_coef`, the directions.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A square base matrix of side `pen@n_coef`, a list of them keyed by
#'   `pen@params`, or an error.
#'
#' @seealso [penalty_d2hessian_beta()] for the generics,
#'   [penalty_dhessian_beta.penalty()] for the methods of the derivative of
#'   first order.
#' @keywords internal
S7::method(penalty_d2hessian_beta, penalty) <- function(pen, beta, theta, v, w,
                                                        ...) {
  stop(sprintf(paste0("'%s' does not supply penalty_d2hessian_beta(), so the\n",
                      "  second movement of its Hessian with the coefficients",
                      " is not available."),
               pen@penalty_name), call. = FALSE)
}

#' @rdname penalty_d2hessian_beta.penalty
#' @name penalty_dhessian_beta_theta.penalty
#' @keywords internal
S7::method(penalty_dhessian_beta_theta, penalty) <- function(pen, beta, theta,
                                                             v, ...) {
  stop(sprintf(paste0("'%s' does not supply penalty_dhessian_beta_theta(), so",
                      "\n  how the movement of its Hessian changes with the",
                      " hyperparameters is not available."),
               pen@penalty_name), call. = FALSE)
}

# the three branches whose Hessian does not depend on the coefficients; the
# methods are documented on the page penalty_d2hessian_beta.penalty, which
# carries their aliases
for (.cls in list(QuadraticPenalty, AdditivePenalty, StructuredPenalty)) {
  S7::method(penalty_d2hessian_beta, .cls) <- function(pen, beta, theta, v, w,
                                                       ...) {
    matrix(0, as.integer(pen@n_coef), as.integer(pen@n_coef))
  }
  S7::method(penalty_dhessian_beta_theta, .cls) <- function(pen, beta, theta,
                                                            v, ...) {
    k <- as.integer(pen@n_coef)
    stats::setNames(lapply(pen@params, function(m) matrix(0, k, k)),
                    pen@params)
  }
}
rm(.cls)

#' @rdname penalty_d2hessian_beta.penalty
#' @name penalty_d2hessian_beta.DistribPenalty
#' @keywords internal
S7::method(penalty_d2hessian_beta, DistribPenalty) <- function(pen, beta,
                                                               theta, v, w,
                                                               ...) {
  reject_kinked(pen, "penalty_d2hessian_beta")
  k <- as.integer(pen@n_coef)
  if (isTRUE(beta_quadratic(pen, theta))) return(matrix(0, k, k))
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    d4 <- dp_parent_or_reject(pen, "fourth", distributions7::distrib_deriv4_y,
                              a, theta)
    h <- dp_contract(d4, list(dp_arg(pen, map_apply(pen, v)),
                              dp_arg(pen, map_apply(pen, w))))
    return(-map_quad_full(pen, dp_blockdiag(pen, h, nrow(a))))
  }
  d4 <- distributions7::distrib_deriv4_y(pen@parent, t, theta) + 0 * t
  -map_quad(pen, d4 * map_apply(pen, v) * map_apply(pen, w))
}

#' @rdname penalty_d2hessian_beta.penalty
#' @name penalty_dhessian_beta_theta.DistribPenalty
#' @keywords internal
S7::method(penalty_dhessian_beta_theta, DistribPenalty) <- function(pen, beta,
                                                                    theta, v,
                                                                    ...) {
  reject_kinked(pen, "penalty_dhessian_beta_theta")
  k <- as.integer(pen@n_coef)
  if (isTRUE(beta_quadratic(pen, theta))) {
    return(stats::setNames(lapply(pen@params, function(m) matrix(0, k, k)),
                           pen@params))
  }
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    c3 <- dp_parent_or_reject(pen, "mixed third",
                              distributions7::distrib_cross3_y, a, theta)
    av <- dp_arg(pen, map_apply(pen, v))
    return(stats::setNames(lapply(pen@params, function(m)
      -map_quad_full(pen, dp_blockdiag(pen, dp_contract(c3[[m]], list(av)),
                                       nrow(a)))), pen@params))
  }
  tv <- map_apply(pen, v)
  c3 <- distributions7::distrib_cross3_y(pen@parent, t, theta)
  stats::setNames(lapply(pen@params, function(m)
    -map_quad(pen, (c3[[m]] + 0 * t) * tv)), pen@params)
}


#' Contracting a Multivariate Parent's Response Tensor Per Block
#'
#' @description
#' `dp_contract()` contracts the trailing response indices of a per-block
#' tensor against one direction each: given an array of dimension
#' \eqn{p \times \dots \times p \times n} and \eqn{m} matrices of \eqn{n}
#' rows and \eqn{p} columns, it returns the \eqn{p \times p \times n} array
#' whose slice \eqn{i} is
#' \eqn{\sum_{c_1 \dots c_m} A_{ab c_1 \dots c_m i}\,v_{1,ic_1}\cdots
#' v_{m,ic_m}}. That is the derivative of the response Hessian of block
#' \eqn{i} along the coordinates of the block in the directions, which
#' [dp_blockdiag()] then places on the diagonal. Arrays of four and five
#' dimensions are handled, that is one or two directions.
#'
#' `dp_parent_or_reject()` calls the generic `fun` of the parent. If that call
#' signals an error of any kind, it signals in its place an error that names
#' the penalty and the derivative, so that a criterion receives either the
#' complete tensor or an error.
#'
#' @param A A numeric array with \eqn{n} as its last dimension.
#' @param dirs A list of \eqn{n \times p} matrices, one per contracted index;
#'   the last remaining index of the array is contracted first, against the
#'   last direction.
#' @param pen A [DistribPenalty()] object.
#' @param what A word naming the derivative, for the message.
#' @param fun The \pkg{distributions7} generic to call.
#' @param a The argument of the parent, from [dp_arg()].
#' @param theta A named list of the free parameters of the parent.
#'
#' @return `dp_contract()` a \eqn{p \times p \times n} numeric array.
#'   `dp_parent_or_reject()` the value returned by `fun`, or an error.
#'
#' @seealso [penalty_dhessian_beta()], [penalty_d2hessian_beta()].
#'
#' @examples
#' A <- array(1, c(2, 2, 2, 3))
#' v <- matrix(c(1, 0, 2, 1, 1, 0), 3, 2)
#' dim(penalties7:::dp_contract(A, list(v)))
#'
#' @keywords internal
dp_contract <- function(A, dirs) {
  d <- dim(A)
  n <- d[length(d)]
  p <- d[1L]
  out <- array(0, c(p, p, n))
  for (i in seq_len(n)) {
    s <- switch(as.character(length(d)),
                "4" = A[, , , i], "5" = A[, , , , i])
    for (k in rev(seq_along(dirs))) {
      dk <- length(dim(s))
      s <- apply(s, seq_len(dk - 1L), function(x) sum(x * dirs[[k]][i, ]))
    }
    out[, , i] <- s
  }
  out
}

#' @rdname dp_contract
#' @keywords internal
dp_parent_or_reject <- function(pen, what, fun, a, theta) {
  tryCatch(fun(pen@parent, a, theta), error = function(e) {
    stop(sprintf(paste0("'%s' has a multivariate parent that is not quadratic,",
                        "\n  and its %s response derivative per block is not",
                        " available."),
                 pen@penalty_name, what), call. = FALSE)
  })
}
