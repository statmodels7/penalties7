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
#' \partial\theta_m}, one matrix per hyperparameter, saying how the penalty's
#' curvature in the coefficients moves as each hyperparameter moves. This is
#' the third-order quantity a marginal likelihood criterion needs and that the
#' value, the gradient and the Hessian cannot supply between them.
#'
#' @details
#' # Why a criterion needs it
#'
#' Estimating hyperparameters by a marginal likelihood means differentiating a
#' Laplace approximation, whose determinant is \eqn{\lvert H + S\rvert} with
#' \eqn{S} the penalty's Hessian in the coefficients. Differentiating
#' \eqn{\log\lvert H+S\rvert} in \eqn{\theta_m} needs
#' \eqn{\partial S/\partial\theta_m}, and nothing below third order gives it.
#'
#' # What each branch answers
#'
#' | branch | \eqn{\partial S/\partial\theta_m} |
#' |---|---|
#' | [quadratic_penalty()] | \eqn{D'PD}, free of both \eqn{\lambda} and \eqn{\beta} |
#' | [additive_penalty()] | the component \eqn{P_k} of each smoothing parameter |
#' | [structured_penalty()] | the matrix parameter's own `param_d1` |
#' | [distrib_penalty()] | \eqn{-D'\mathrm{diag}(\partial^3\ell/\partial y^2\partial\theta_m)D}, from [distributions7::distrib_cross2_y()] |
#' | [scad_penalty()], [mcp_penalty()] | rejects |
#'
#' Every one of the first four is closed form. The base class rejects, naming
#' what is missing, so a penalty written later either supplies the three
#' methods and works in a marginal criterion, or says plainly that it cannot.
#'
#' A separable penalty with a **kink** rejects as well, whatever its parent:
#' the third derivative does not exist there, and the mode a marginal criterion
#' expands around is exactly where a kinked penalty puts coefficients.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Read only by the
#'   separable branch, whose parent's derivatives depend on the argument; the
#'   quadratic, additive and structured branches ignore it.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A named list of one matrix per hyperparameter, each of side
#'   `pen@n_coef`, named by `pen@params`. The quadratic branch returns the
#'   stored matrix in whatever storage it keeps, so a penalty built with
#'   `blocks > 1` gives a `dgCMatrix`.
#'
#' @seealso [penalty_hessian()] for the quantity differentiated,
#'   [penalty_d2hessian()] and [penalty_dcross()] for the other two a marginal
#'   criterion asks for, [beta_quadratic()] for the predicate that says whether
#'   a third \eqn{\beta}-derivative is needed at all.
#'
#' @examples
#' # A quadratic penalty's Hessian is lambda D'PD, so the derivative is D'PD
#' # and does not move with lambda.
#' pen <- quadratic_penalty(crossprod(diff(diag(3))))
#' b <- c(1, -0.5, 0.3)
#' penalty_dhessian(pen, b, list(lambda = 2))$lambda
#' identical(penalty_dhessian(pen, b, list(lambda = 2)),
#'           penalty_dhessian(pen, b, list(lambda = 99)))
#'
#' # An additive penalty answers with its components, one per parameter.
#' add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
#' names(penalty_dhessian(add, b, list(lambda1 = 2, lambda2 = 0.5)))
#'
#' # A penalty with a kink has no such quantity and says so.
#' try(penalty_dhessian(lasso_penalty(n_coef = 3), b, list(lambda = 1)))
#'
#' @export
penalty_dhessian <- S7::new_generic("penalty_dhessian", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    S7::S7_dispatch()
  })


#' The Second Derivative of the Coefficient Hessian in the Hyperparameters
#'
#' @description
#' Returns \eqn{\partial^2 S/\partial\theta_m\partial\theta_l =
#' \partial^4\rho/\partial\beta^2\partial\theta_m\partial\theta_l}, one matrix
#' per unordered pair. A marginal criterion differentiated twice in the
#' hyperparameters needs it, and so does an exact outer Hessian.
#'
#' @details
#' It is **exactly zero** for every penalty whose Hessian is linear in its
#' hyperparameters, which is the quadratic branch, where
#' \eqn{S = \lambda D'PD}, and the additive branch, where
#' \eqn{S = \sum_k \lambda_k P_k}. The structured branch reads the matrix
#' parameter's `param_d2`, and the separable branch the second
#' \eqn{\theta}-derivative of the parent's response curvature, through
#' [distributions7::distrib_hess_y_hess()].
#'
#' The keys are [penalty_hess_theta()]'s: diagonals first in `pen@params`
#' order, then the upper off-diagonal pairs, each joined by an underscore. A
#' consumer looking a pair up therefore need not know which order it was
#' written in.
#'
#' A penalty with a kink rejects, and so does the base class.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Read only by the
#'   separable branch.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same.
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
#' # And for an additive one, keyed diagonals first over three pairs.
#' add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
#' names(penalty_d2hessian(add, b, list(lambda1 = 2, lambda2 = 0.5)))
#'
#' # Not zero for a structured penalty, whose matrix is not linear in any
#' # one free value.
#' s <- structured_penalty(parameters7::log_cholesky(2))
#' th <- list(log_L1 = 0.1, log_L2 = -0.1, L2.1 = 0.3)
#' penalty_d2hessian(s, c(1, -0.5), th)$log_L1_log_L1
#'
#' @export
penalty_d2hessian <- S7::new_generic("penalty_d2hessian", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    S7::S7_dispatch()
  })


#' The Derivative of the Mixed Block in the Hyperparameters
#'
#' @description
#' Returns \eqn{\partial^3\rho/\partial\beta\,\partial\theta_m\partial\theta_l},
#' one coefficient vector per unordered pair. Where [penalty_cross()] says how
#' the coefficient gradient moves with one hyperparameter, this says how that
#' movement itself moves with a second.
#'
#' @details
#' A marginal criterion is evaluated at the penalized mode, and the mode moves
#' with the hyperparameters. Differentiating the criterion twice therefore
#' needs the derivative of the mode's own derivative, and this is the penalty's
#' contribution to it.
#'
#' It is **exactly zero** wherever the penalty is quadratic in the coefficients
#' with a Hessian linear in the hyperparameters, which covers the quadratic and
#' additive branches. The structured branch answers
#' \eqn{(\partial^2\Omega/\partial\theta_m\partial\theta_l)\beta}, reusing
#' [penalty_d2hessian()], and the separable branch reads the parent's
#' [distributions7::distrib_grad_y_hess()].
#'
#' The keys are [penalty_hess_theta()]'s, diagonals first.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Read by the structured
#'   and separable branches; the quadratic and additive ones answer zero
#'   whatever it is.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A named list of one numeric vector of length `pen@n_coef` per
#'   unordered pair, keyed diagonals first.
#'
#' @seealso [penalty_cross()] for the first order,
#'   [penalty_dhessian()] and [penalty_d2hessian()] for the other two
#'   quantities a marginal criterion asks for.
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
    S7::S7_dispatch()
  })


#' Is a Penalty Quadratic in the Coefficients?
#'
#' @description
#' `TRUE` when \eqn{\partial^3\rho/\partial\beta^3} is exactly zero, so a
#' consumer can skip that third derivative altogether. Asked of the penalty
#' rather than measured by a consumer, which would be guessing at a property
#' the penalty knows.
#'
#' @details
#' # Three independent properties
#'
#' This is not [is_quadratic()], which is about the whole construction. The
#' three questions a consumer may ask are independent:
#'
#' | penalty | `is_quadratic()` | `beta_quadratic()` | \eqn{S} linear in \eqn{\theta} |
#' |---|---|---|---|
#' | [quadratic_penalty()] | `TRUE` | `TRUE` | yes |
#' | [structured_penalty()] | `TRUE` | `TRUE` | no |
#' | [additive_penalty()] | `FALSE` | `TRUE` | yes |
#' | a Gaussian [distrib_penalty()] | `FALSE` | `TRUE` | no |
#' | [heavy_penalty()] | `FALSE` | `FALSE` | no |
#'
#' # How the separable branch answers
#'
#' By probing, not by differentiating three times. The log-density is quadratic
#' in the response exactly when its second derivative there does not depend on
#' it, so [distributions7::distrib_hess_y()] is read at four points and
#' compared for constancy at a tolerance of `1e-12`. The question is put to the
#' second derivative because that one is analytic for almost every family,
#' where the third is often a finite difference whose noise no threshold
#' separates from a true zero.
#'
#' # What `TRUE` does not promise
#'
#' It does not say the marginal derivatives are available. A lasso answers
#' `TRUE`, its log-density being linear in the response away from the kink, and
#' [penalty_dhessian()] still rejects for it: the third derivative is zero
#' where it exists and undefined at the kink, which is where a selecting fit
#' puts its coefficients.
#'
#' @param pen A [penalty()] object.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same. Read by the separable branch, to evaluate the
#'   parent; the other branches answer a constant.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A single logical. `FALSE` from the base class, so a penalty that
#'   says nothing is treated as needing the third derivative.
#'
#' @seealso [is_quadratic()] for the different question about the
#'   construction, [penalty_dhessian()] for what a consumer skips when this is
#'   `TRUE`, [penalty_hessian()] for the derivative that is then constant.
#'
#' @examples
#' # The two branches that answer TRUE to both predicates.
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
  function(pen, theta, ...) S7::S7_dispatch())


# --- the base class rejects -------------------------------------------------

#' @title What the Base Class Answers to the Marginal Generics
#' @name penalty_dhessian.penalty
#'
#' @description
#' One page for the base class's answers to the four generics a marginal
#' criterion asks beyond the second order. `penalty_dhessian()`,
#' `penalty_d2hessian()` and `penalty_dcross()` each reject, naming the penalty
#' and the generic; `beta_quadratic()` answers `FALSE`.
#'
#' @details
#' A marginal criterion is a Laplace approximation and asks for derivatives
#' beyond the second. A penalty that does not supply them cannot be estimated
#' by one, and reporting that is better than a criterion assembled from a
#' quantity nobody wrote. `penalty_dhessian()`'s message says so in full; the
#' other two name the generic alone, a caller who has reached the second having
#' already read the first.
#'
#' `beta_quadratic()`'s `FALSE` is the conservative default: a penalty that
#' says nothing about its third \eqn{\beta}-derivative is treated as having
#' one.
#'
#' @param pen A [penalty()] object of a class registering none of these.
#' @param beta A numeric vector of coefficients. Unused.
#' @param theta A named list of hyperparameter values. Unused.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_dhessian()`, `penalty_d2hessian()` and `penalty_dcross()`
#'   signal an error. `beta_quadratic()` returns the single logical `FALSE`.
#'
#' @seealso [penalty_dhessian()] and its siblings for the generics,
#'   [penalty_dhessian.QuadraticPenalty()] for a branch that answers.
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
#' One page for the branch's answers to the four generics a marginal criterion
#' asks beyond the second order. \eqn{S = \lambda D'PD} is linear in
#' \eqn{\lambda} and free of the coefficients, so `penalty_dhessian()` returns
#' \eqn{D'PD}, `penalty_d2hessian()` and `penalty_dcross()` return zeros, and
#' `beta_quadratic()` returns `TRUE`.
#'
#' @details
#' `penalty_dhessian()` returns the stored \eqn{D'PD} **in whatever storage the
#' penalty keeps it**, `unclass()`ed to strip the attributes a base matrix may
#' carry and leaving an S4 matrix untouched. A penalty built with `blocks > 1`
#' therefore answers with a `dgCMatrix`, and a consumer that needs a base
#' matrix coerces where the two meet.
#'
#' The zeros are exact. The Hessian is \eqn{\lambda \cdot} a constant, so its
#' second derivative in \eqn{\lambda} vanishes identically, and so does the
#' third derivative of the value in \eqn{\beta} and two hyperparameters.
#'
#' @param pen A [QuadraticPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Unused: no quantity
#'   here depends on the coefficients.
#' @param theta A named list holding `lambda`. Unused for the same reason.
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
#' # The derivative is D'PD, so it does not move with lambda at all.
#' max(abs(penalty_dhessian(pen, b, list(lambda = 2))$lambda -
#'           penalty_dhessian(pen, b, list(lambda = 99))$lambda))
#'
#' # And everything above first order is exactly zero.
#' max(abs(unlist(penalty_d2hessian(pen, b, list(lambda = 2)))))
#' max(abs(unlist(penalty_dcross(pen, b, list(lambda = 2)))))
#' beta_quadratic(pen, list(lambda = 2))
#'
#' @seealso [penalty_dhessian()] and its siblings for the generics,
#'   [quadratic_penalty()] for the branch, [zero_pairs()] for the zeros.
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
#' One page for the branch's answers to the four generics a marginal criterion
#' asks beyond the second order. \eqn{S = \sum_k \lambda_k P_k} is linear in
#' the smoothing parameters and free of the coefficients, so
#' `penalty_dhessian()` returns the components themselves,
#' `penalty_d2hessian()` and `penalty_dcross()` return zeros, and
#' `beta_quadratic()` returns `TRUE`.
#'
#' @details
#' Every entry is exact and none is computed at call time: the components were
#' fixed at construction, already carried through the map, so the first
#' derivative is a lookup and the higher ones are zeros of the right shape.
#'
#' This is the one branch whose marginal derivatives are reachable while
#' [is_quadratic()] answers `FALSE` for it. A consumer routing on that
#' predicate rather than on these methods will not find them; see
#' [additive_penalty()].
#'
#' @param pen An [AdditivePenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Unused.
#' @param theta A named list of `lambda1`, `lambda2`, ... Unused.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_dhessian()` a list of one matrix per component, named by
#'   `pen@params`.
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
#'   [additive_penalty()] for the branch and for what `is_quadratic()` costs
#'   it.
#' @keywords internal
S7::method(penalty_dhessian, AdditivePenalty) <- function(pen, beta, theta,
                                                          ...) {
  stats::setNames(lapply(pen@mats, function(P) unclass(P)), pen@params)
}

#' @rdname penalty_dhessian.AdditivePenalty
#' @name penalty_d2hessian.AdditivePenalty
#' @keywords internal
S7::method(penalty_d2hessian, AdditivePenalty) <- function(pen, beta, theta,
                                                           ...) {
  zero_pairs(pen, matrix(0, as.integer(pen@n_coef), as.integer(pen@n_coef)))
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
#' One page for the branch's answers to the four generics a marginal criterion
#' asks beyond the second order. \eqn{S = \Omega(\theta)} is the matrix
#' parameter itself, so `penalty_dhessian()` and `penalty_d2hessian()` are the
#' structure's own `param_d1` and `param_d2`, `penalty_dcross()` is the second
#' of those applied to the coefficients, and `beta_quadratic()` is `TRUE`.
#'
#' @details
#' Nothing is differentiated here: \pkg{parameters7} supplies both derivative
#' arrays exactly for every structure it ships, and this branch unwraps and
#' re-keys them.
#'
#' Since the value is quadratic in the coefficients, the mixed third derivative
#' is \eqn{(\partial^2\Omega/\partial\theta_m\partial\theta_l)\beta}, so
#' `penalty_dcross()` calls `penalty_d2hessian()` and multiplies.
#'
#' # The keying
#'
#' `param_d2` is a **flat** list whose keys join the free names with a
#' separator of its own, so an entry is located through
#' [parameters7::param_tuple_indices()], which enumerates the tuples in exactly
#' the order the components are in. Nothing here depends on the spelling of a
#' key, which matters because a free name may itself contain the separator. A
#' pair with no component raises rather than returning something of the wrong
#' shape.
#'
#' @param pen A [StructuredPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Read by
#'   `penalty_dcross()`; the other three ignore it.
#' @param theta A named list of the structure's free values.
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
#'   [parameters7::param_d1()] for the arrays it reads.
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
#' One page for the branch's answers to the four generics a marginal criterion
#' asks beyond the second order. With
#' \eqn{\rho = -\sum_j \log f((D\beta)_j;\theta)} the Hessian is
#' \eqn{-D'\mathrm{diag}(\ell^{(yy)})D}, so its \eqn{\theta}-derivatives are
#' the parent's own higher components carried through the same map.
#'
#' @details
#' # Where each comes from
#'
#' | generic | parent's component |
#' |---|---|
#' | `penalty_dhessian()` | [distributions7::distrib_cross2_y()] |
#' | `penalty_d2hessian()` | [distributions7::distrib_hess_y_hess()] |
#' | `penalty_dcross()` | [distributions7::distrib_grad_y_hess()] |
#'
#' Nothing is differentiated here. A parent with closed forms for those, the
#' gaussian among them and so every Gaussian random effect, makes this branch
#' exact; a parent without them inherits that package's documented fallback,
#' one central difference of its analytic first-order component.
#'
#' A multivariate parent is read blockwise: each component is assembled into a
#' block-diagonal matrix by [dp_blockdiag()] and carried through the map by
#' [map_quad_full()], where a univariate one goes through [map_quad()].
#'
#' # What rejects
#'
#' A penalty whose parent declares a kink, which is the lasso and the elastic
#' net. The third derivative does not exist at the kink, and the mode a
#' marginal criterion expands around is exactly where a kinked penalty puts
#' coefficients, so all three derivative generics reject with the penalty
#' named. `beta_quadratic()` does not: it answers `TRUE` for the lasso, whose
#' log-density is linear in the response away from the kink.
#'
#' # How `beta_quadratic()` decides
#'
#' By probing the parent's [distributions7::distrib_hess_y()] at four points
#' and asking whether it is constant to `1e-12`. The log-density is quadratic
#' in the response exactly when its second derivative there does not depend on
#' it, and that second derivative is analytic for almost every family, where
#' the third is often a difference whose noise no threshold separates from a
#' true zero.
#'
#' @param pen A [DistribPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`. Read by all three
#'   derivative methods, the parent's components depending on the argument.
#' @param theta A named list of the parent's free parameters.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_dhessian()` a list of one matrix of side `pen@n_coef` per
#'   hyperparameter.
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
#' # derivative moves with them.
#' h <- heavy_penalty(n_coef = 3)
#' beta_quadratic(h, list(sigma = 1, nu = 4))
#' diag(penalty_dhessian(h, b, list(sigma = 1, nu = 4))$nu)
#'
#' # A kinked parent has no third derivative and says so.
#' try(penalty_dhessian(lasso_penalty(n_coef = 3), b, list(lambda = 1)))
#'
#' @seealso [penalty_dhessian()] and its siblings for the generics,
#'   [distrib_penalty()] for the branch, [reject_kinked()] for the refusal,
#'   [carry_pairs()] for the re-keying.
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
  t <- c(-1.73, -0.29, 0.61, 2.04)
  h <- tryCatch(distributions7::distrib_hess_y(pen@parent, t,
                                               align_ptheta(pen, theta)),
                error = function(e) NULL)
  !is.null(h) && length(h) > 1L &&
    isTRUE(all.equal(as.numeric(h), rep(as.numeric(h)[1L], length(h)),
                     tolerance = 1e-12))
}


#' Carry a Parent's Paired Components Into Coefficient Space
#'
#' @description
#' Re-keys a \pkg{distributions7} component keyed by parameter pair into this
#' package's own keys, and places each through the penalty's map.
#'
#' @details
#' The two enumerations of pairs are built the same way from the same names, so
#' a key from one is a key of the other; it is looked up by name in both
#' orders rather than by position, since a hyperparameter whose own name
#' contains the separator would not survive being taken apart.
#'
#' @param pen A [penalty()] object.
#' @param comp The parent's components, keyed by parameter pair.
#' @param carry A function placing one component into coefficient space.
#'
#' @return A named list keyed by hyperparameter pair.
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
#' key, which is the answer of a penalty whose Hessian is linear in its
#' hyperparameters.
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
#' Signals an error naming the penalty and the generic when the penalty carries
#' a kink, and returns invisibly otherwise. Called at the head of each of the
#' separable branch's three marginal derivative methods.
#'
#' @details
#' A marginal criterion is a Laplace expansion at the penalized mode, and a
#' kinked penalty puts coefficients exactly at the kink, where the third
#' derivative does not exist. Rejecting there is better than returning the
#' one-sided value, which would give a criterion no one could interpret.
#'
#' The kinks are read straight off `pen@kinks` rather than through
#' [penalty_kinks()], which takes a `theta` this is asked before there is one
#' to pass. An object without that property, which no shipped branch is, is
#' treated as having none.
#'
#' @param pen A [penalty()] object.
#' @param what The generic's name, a single string, quoted back in the message.
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
