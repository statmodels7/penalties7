#' @include penalty_class.R
NULL

#' @title Value of a Penalty
#'
#' @description
#' The scalar \eqn{\rho(D\beta;\theta)}, with the normalizing constant
#' included whenever the penalty is proper, so that the value is exactly the
#' negative log-density of the prior.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return A single number.
#'
#' @examples
#' pen <- quadratic_penalty(diag(3))
#' penalty_value(pen, c(1, 0, -1), list(lambda = 2))
#'
#' @seealso [penalty_gradient()], [penalty_hessian()], [penalty_grad_theta()], [penalty_cross()], [penalty_kinks()]
#' @export
penalty_value <- S7::new_generic("penalty_value", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    S7::S7_dispatch()
  })

#' @title Coefficient Derivatives of a Penalty
#'
#' @description
#' `penalty_gradient` returns \eqn{\partial\rho/\partial\beta} and
#' `penalty_hessian` returns \eqn{\partial^2\rho/\partial\beta^2}, both
#' exact.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return `penalty_gradient` a numeric vector of length `q`;
#'   `penalty_hessian` a `q x q` symmetric matrix.
#'
#' @examples
#' pen <- quadratic_penalty(diag(2))
#' penalty_gradient(pen, c(1, -1), list(lambda = 3))
#' penalty_hessian(pen, c(1, -1), list(lambda = 3))
#'
#' @seealso [penalty_value()], [penalty_grad_theta()], [penalty_kinks()]
#' @export
penalty_gradient <- S7::new_generic("penalty_gradient", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    S7::S7_dispatch()
  })

#' @rdname penalty_gradient
#' @export
penalty_hessian <- S7::new_generic("penalty_hessian", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    S7::S7_dispatch()
  })

#' @title Hyperparameter Derivatives of a Penalty
#'
#' @description
#' `penalty_grad_theta` returns \eqn{\partial\rho/\partial\theta} as a
#' named list, `penalty_hess_theta` the second derivatives keyed with
#' diagonals first, and `penalty_cross` the mixed block
#' \eqn{\partial^2\rho/\partial\beta\,\partial\theta_k}, one coefficient
#' vector per hyperparameter -- the block a joint estimation of coefficients
#' and hyperparameters needs.
#'
#' @param pen A [penalty()] object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameter values.
#' @param scale Either `"parameter"` (default) or `"link"`; on the
#'   link scale the derivatives are with respect to the unconstrained
#'   values, carried by the chain rule in the generic body, so methods
#'   always return the parameter scale.
#' @param ... Passed to methods.
#'
#' @return A named list: one number per hyperparameter for the gradient, one
#'   number per pair for the Hessian, one numeric vector of length `q`
#'   per hyperparameter for the mixed block.
#'
#' @examples
#' pen <- quadratic_penalty(diag(2))
#' penalty_grad_theta(pen, c(1, -1), list(lambda = 3))
#' penalty_cross(pen, c(1, -1), list(lambda = 3))
#'
#' @seealso [penalty_value()], [penalty_gradient()], [penalty_kinks()]
#' @export
penalty_grad_theta <- S7::new_generic("penalty_grad_theta", "pen",
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    scale <- match.arg(scale)
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    out <- S7::S7_dispatch()
    if (scale == "link") out <- ptheta_to_link(pen, theta, g = out)
    out
  })

#' @rdname penalty_grad_theta
#' @export
penalty_hess_theta <- S7::new_generic("penalty_hess_theta", "pen",
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    scale <- match.arg(scale)
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    out <- S7::S7_dispatch()
    if (scale == "link") {
      g <- penalty_grad_theta(pen, beta, theta)
      out <- ptheta_to_link(pen, theta, g = g, H = out)
    }
    out
  })

#' @rdname penalty_grad_theta
#' @export
penalty_cross <- S7::new_generic("penalty_cross", "pen",
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    scale <- match.arg(scale)
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    out <- S7::S7_dispatch()
    if (scale == "link") out <- ptheta_to_link(pen, theta, cross = out)
    out
  })

#' @title The Non-Differentiable Points of a Penalty
#'
#' @description
#' The values of \eqn{t = D\beta} at which \eqn{\rho} is not differentiable
#' in its argument: empty for the smooth penalties, `0` for the lasso,
#' SCAD and MCP. [check_penalty()] places its grids away from
#' them, and a solver may consult them.
#'
#' @param pen A [penalty()] object.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return A numeric vector, possibly empty.
#'
#' @examples
#' penalty_kinks(quadratic_penalty(diag(2)), list(lambda = 1))
#' penalty_kinks(lasso_penalty(), list(lambda = 1))
#'
#' @seealso [penalty_value()], [penalty_gradient()], [penalty_hessian()], [penalty_grad_theta()], [penalty_cross()]
#' @export
penalty_kinks <- S7::new_generic("penalty_kinks", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

#' @title Is a Penalty a Proper Prior?
#'
#' @description
#' `TRUE` when \eqn{\exp(-\rho)} integrates to one over the penalized
#' coordinates, so that the value is exactly a negative log-density;
#' `FALSE` for the improper ones (a rank-deficient quadratic, SCAD,
#' MCP), whose value is the bare \eqn{\rho}.
#'
#' @param pen A [penalty()] object.
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @examples
#' is_proper(quadratic_penalty(diag(2)))
#' is_proper(scad_penalty())
#'
#' @seealso [is_quadratic()], [has_prox()], [penalty_matrix()], [penalty_rank()], [penalty_null_basis()], [penalty_logpdet()]
#' @export
is_proper <- S7::new_generic("is_proper", "pen")

#' @title Is a Penalty Quadratic?
#'
#' @description
#' `TRUE` only for [quadratic_penalty()], whose matrix, rank,
#' null basis and log pseudo-determinant the marginal-likelihood generics
#' expose.
#'
#' @details
#' A quadratic penalty carries a fixed matrix \eqn{P} and one smoothing
#' parameter \eqn{\lambda}, and is the negative log-density of the improper
#' Gaussian prior with precision \eqn{\lambda P} on \eqn{D\beta}:
#'
#' \deqn{\rho(\beta; \lambda)
#'   = \tfrac{\lambda}{2} (D\beta)^\top P (D\beta)
#'   - \tfrac{1}{2}\log^{+}\lvert \lambda P \rvert
#'   + \tfrac{r}{2}\log(2\pi),
#'   \qquad \log^{+}\lvert \lambda P \rvert
#'     = r \log \lambda + \log^{+}\lvert P \rvert,}
#'
#' with \eqn{\log^{+}} the log pseudo-determinant and
#' \eqn{r = \operatorname{rank}(P)}. Those are the quantities
#' [penalty_matrix()], [penalty_rank()],
#' [penalty_null_basis()] and [penalty_logpdet()]
#' report and a REML or marginal-likelihood criterion needs; a penalty for
#' which this is `FALSE` has no such matrix and those generics reject.
#'
#' @param pen A [penalty()] object.
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @examples
#' is_quadratic(quadratic_penalty(diag(2)))
#' is_quadratic(ridge_penalty())
#'
#' @seealso [is_proper()], [has_prox()], [penalty_matrix()], [penalty_rank()], [penalty_null_basis()], [penalty_logpdet()]
#' @export
is_quadratic <- S7::new_generic("is_quadratic", "pen")

S7::method(is_quadratic, penalty) <- function(pen, ...) FALSE

#' @title The Pieces a Marginal Criterion Consumes
#'
#' @description
#' For a quadratic penalty: `penalty_matrix` returns
#' \eqn{\lambda D'PD}, `penalty_rank` its rank (fixed at construction),
#' `penalty_null_basis` the exact null basis for the model layer to
#' intersect across terms, and `penalty_logpdet` the log
#' pseudo-determinant \eqn{r\log\lambda + \log\mathrm{pdet}(P)} with its
#' first two theta derivatives. Every other penalty rejects: a marginal
#' criterion for a non-Gaussian prior is not a determinant, and pretending
#' otherwise would produce numbers silently.
#'
#' @param pen A [penalty()] object.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return `penalty_matrix` a `q x q` matrix; `penalty_rank`
#'   an integer; `penalty_null_basis` a matrix with `q` rows (zero
#'   columns when the penalty is full rank); `penalty_logpdet` a list
#'   with elements `value`, `grad` and `hess`.
#'
#' @examples
#' pen <- quadratic_penalty(crossprod(diff(diag(4))))
#' penalty_rank(pen)
#' penalty_logpdet(pen, list(lambda = 2))$value
#'
#' @seealso [is_quadratic()], [quadratic_penalty()], [additive_penalty()]
#' @export
penalty_matrix <- S7::new_generic("penalty_matrix", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

#' @rdname penalty_matrix
#' @export
penalty_rank <- S7::new_generic("penalty_rank", "pen")

#' @rdname penalty_matrix
#' @export
penalty_null_basis <- S7::new_generic("penalty_null_basis", "pen")

#' @rdname penalty_matrix
#' @export
penalty_logpdet <- S7::new_generic("penalty_logpdet", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

S7::method(penalty_matrix, penalty) <- function(pen, theta, ...) {
  stop("Only a quadratic penalty exposes its matrix; see is_quadratic().",
       call. = FALSE)
}
S7::method(penalty_rank, penalty) <- function(pen, ...) {
  stop("Only a quadratic penalty has a rank; see is_quadratic().",
       call. = FALSE)
}
S7::method(penalty_null_basis, penalty) <- function(pen, ...) {
  stop("Only a quadratic penalty has a null basis; see is_quadratic().",
       call. = FALSE)
}
S7::method(penalty_logpdet, penalty) <- function(pen, theta, ...) {
  stop("Only a quadratic penalty has a log pseudo-determinant; see is_quadratic().",
       call. = FALSE)
}


#' @title What a Penalty's Hyperparameters Are About
#'
#' @description
#' The quantities a reader reads, where the hyperparameters are coordinates of
#' a chart rather than the quantities themselves, with the Jacobian from those
#' coordinates and the scale each one's interval belongs on.
#'
#' @details
#' The case this exists for is a penalty whose prior is a multivariate family:
#' its hyperparameters are the free values of a matrix parameter -- the
#' logarithms of the diagonal of a Cholesky factor and the entries below it --
#' and nobody reads those. What the prior is about is the standard deviations
#' and the correlations of the effects it describes, and
#' [distributions7::mv_derived()] declares them, so this is the same
#' distinction [parameters7::param_readable()] makes for a matrix
#' parameter and `term_readable` for a fitted term.
#'
#' The base method returns `NULL`, which says that the hyperparameters
#' ARE the quantities and a consumer should report them as they stand. That is
#' the honest answer for every other branch: a smoothing parameter, a rate, a
#' shape are each read on their own scale already.
#'
#' @param pen A [penalty()] object.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return `NULL`, or a list with `value`, `jacobian`,
#'   `transform` and `block`, as
#'   [distributions7::mv_derived()] returns them.
#'
#' @examples
#' pen <- distrib_penalty(
#'   distributions7::fixed(distributions7::mvgaussian_distrib(2),
#'                         mu1 = 0, mu2 = 0), n_coef = 6)
#' penalty_readable(pen, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
#'                            sigma_L2.1 = 0.5))$value
#'
#' # a smoothing parameter is already the quantity it names
#' penalty_readable(quadratic_penalty(diag(2)), list(lambda = 1))
#'
#' @seealso [penalty_value()], [distributions7::mv_derived()]
#' @export
penalty_readable <- S7::new_generic("penalty_readable", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

S7::method(penalty_readable, penalty) <- function(pen, theta, ...) NULL
