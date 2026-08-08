#' @include penalty_class.R
NULL

#' @title Value of a Penalty
#'
#' @description
#' The scalar \eqn{\rho(D\beta;\theta)}, with the normalizing constant
#' included whenever the penalty is proper, so that the value is exactly the
#' negative log-density of the prior.
#'
#' @param pen A \code{\link{penalty}} object.
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
#' \code{penalty_gradient} returns \eqn{\partial\rho/\partial\beta} and
#' \code{penalty_hessian} returns \eqn{\partial^2\rho/\partial\beta^2}, both
#' exact.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return \code{penalty_gradient} a numeric vector of length \code{q};
#'   \code{penalty_hessian} a \code{q x q} symmetric matrix.
#'
#' @examples
#' pen <- quadratic_penalty(diag(2))
#' penalty_gradient(pen, c(1, -1), list(lambda = 3))
#' penalty_hessian(pen, c(1, -1), list(lambda = 3))
#'
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
#' \code{penalty_grad_theta} returns \eqn{\partial\rho/\partial\theta} as a
#' named list, \code{penalty_hess_theta} the second derivatives keyed with
#' diagonals first, and \code{penalty_cross} the mixed block
#' \eqn{\partial^2\rho/\partial\beta\,\partial\theta_k}, one coefficient
#' vector per hyperparameter -- the block a joint estimation of coefficients
#' and hyperparameters needs.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameter values.
#' @param scale Either \code{"parameter"} (default) or \code{"link"}; on the
#'   link scale the derivatives are with respect to the unconstrained
#'   values, carried by the chain rule in the generic body, so methods
#'   always return the parameter scale.
#' @param ... Passed to methods.
#'
#' @return A named list: one number per hyperparameter for the gradient, one
#'   number per pair for the Hessian, one numeric vector of length \code{q}
#'   per hyperparameter for the mixed block.
#'
#' @examples
#' pen <- quadratic_penalty(diag(2))
#' penalty_grad_theta(pen, c(1, -1), list(lambda = 3))
#' penalty_cross(pen, c(1, -1), list(lambda = 3))
#'
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
#' in its argument: empty for the smooth penalties, \code{0} for the lasso,
#' SCAD and MCP. \code{\link{check_penalty}} places its grids away from
#' them, and a solver may consult them.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return A numeric vector, possibly empty.
#'
#' @examples
#' penalty_kinks(quadratic_penalty(diag(2)), list(lambda = 1))
#' penalty_kinks(lasso_penalty(), list(lambda = 1))
#'
#' @export
penalty_kinks <- S7::new_generic("penalty_kinks", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

#' @title Is a Penalty a Proper Prior?
#'
#' @description
#' \code{TRUE} when \eqn{\exp(-\rho)} integrates to one over the penalized
#' coordinates, so that the value is exactly a negative log-density;
#' \code{FALSE} for the improper ones (a rank-deficient quadratic, SCAD,
#' MCP), whose value is the bare \eqn{\rho}.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @examples
#' is_proper(quadratic_penalty(diag(2)))
#' is_proper(scad_penalty())
#'
#' @export
is_proper <- S7::new_generic("is_proper", "pen")

#' @title Is a Penalty Quadratic?
#'
#' @description
#' \code{TRUE} only for \code{\link{quadratic_penalty}}, whose matrix, rank,
#' null basis and log pseudo-determinant the marginal-likelihood generics
#' expose.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @examples
#' is_quadratic(quadratic_penalty(diag(2)))
#' is_quadratic(ridge_penalty())
#'
#' @export
is_quadratic <- S7::new_generic("is_quadratic", "pen")

S7::method(is_quadratic, penalty) <- function(pen, ...) FALSE

#' @title The Pieces a Marginal Criterion Consumes
#'
#' @description
#' For a quadratic penalty: \code{penalty_matrix} returns
#' \eqn{\lambda D'PD}, \code{penalty_rank} its rank (fixed at construction),
#' \code{penalty_null_basis} the exact null basis for the model layer to
#' intersect across terms, and \code{penalty_logpdet} the log
#' pseudo-determinant \eqn{r\log\lambda + \log\mathrm{pdet}(P)} with its
#' first two theta derivatives. Every other penalty rejects: a marginal
#' criterion for a non-Gaussian prior is not a determinant, and pretending
#' otherwise would produce numbers silently.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param theta A named list of hyperparameter values.
#' @param ... Passed to methods.
#'
#' @return \code{penalty_matrix} a \code{q x q} matrix; \code{penalty_rank}
#'   an integer; \code{penalty_null_basis} a matrix with \code{q} rows (zero
#'   columns when the penalty is full rank); \code{penalty_logpdet} a list
#'   with elements \code{value}, \code{grad} and \code{hess}.
#'
#' @examples
#' pen <- quadratic_penalty(crossprod(diff(diag(4))))
#' penalty_rank(pen)
#' penalty_logpdet(pen, list(lambda = 2))$value
#'
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
