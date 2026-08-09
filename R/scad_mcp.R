#' @include generics.R
NULL

#' @title S7 Classes for the Derivative-Defined Penalties
#'
#' @description
#' The classes \code{\link{scad_penalty}} and \code{\link{mcp_penalty}}
#' instantiate: families the literature defines by \eqn{\rho'} and whose
#' \eqn{\rho} follows by the closed piecewise antiderivative anchored at
#' \eqn{\rho(0) = 0}. Both are improper: \eqn{\rho} is bounded, so
#' \eqn{\exp(-\rho)} does not integrate.
#'
#' @inheritParams penalty
#'
#' @return An object of class \code{ScadPenalty} or \code{McpPenalty}.
#'
#' @seealso \code{\link{scad_penalty}}, \code{\link{mcp_penalty}}
#' @examples
#' S7::S7_inherits(scad_penalty(), ScadPenalty)
#' @keywords internal
#' @export
ScadPenalty <- S7::new_class(name = "ScadPenalty", parent = penalty)

#' @rdname ScadPenalty
#' @keywords internal
#' @export
McpPenalty <- S7::new_class(name = "McpPenalty", parent = penalty)

#' @title Construct the SCAD and MCP Penalties
#'
#' @description
#' SCAD (Fan and Li, 2001) is defined by its derivative on \eqn{t \ge 0},
#' \deqn{\rho'(t) = \lambda \;\; (t \le \lambda), \qquad
#'   \rho'(t) = \dfrac{a\lambda - t}{a - 1} \;\; (\lambda < t \le a\lambda),
#'   \qquad \rho'(t) = 0 \;\; (t > a\lambda),}
#' extended evenly, with \eqn{a > 2}; MCP (Zhang, 2010) by
#' \eqn{\rho'(t) = (\lambda - t/\gamma)_+} with \eqn{\gamma > 1}. Every
#' quantity -- the value, the derivatives in \eqn{\beta} and in the
#' hyperparameters, and the mixed block -- is a closed piecewise form. Both
#' report kinks at zero, where the second derivative additionally jumps at
#' the region boundaries; \code{\link{check_penalty}} keeps its grids away
#' from all of them by asking the object.
#'
#' @param map The matrix \eqn{D}, or \code{NULL} (default) for the identity.
#' @param n_coef The number of coefficients when \code{map} is \code{NULL}.
#' @param link_lambda The link carrying \eqn{\lambda}.
#' @param link_a,link_gamma The link carrying the shoulder parameter, lower
#'   bounded at 2 (SCAD) or 1 (MCP).
#'
#' @return An object of class \code{ScadPenalty} or \code{McpPenalty}.
#'
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. \emph{JASA} 96, 1348-1360.
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. \emph{Annals of Statistics} 38, 894-942.
#'
#' @examples
#' pen <- scad_penalty(n_coef = 3)
#' penalty_value(pen, c(0.5, 2, 5), list(lambda = 1, a = 3.7))
#' is_proper(pen)
#'
#' @export
scad_penalty <- function(map = NULL, n_coef = 1L,
                         link_lambda = linkfunctions7::log_link(),
                         link_a = linkfunctions7::bounded_link(lwr = 2)) {
  q <- if (is.null(map)) as.integer(n_coef) else ncol(map <- as.matrix(map))
  ScadPenalty(
    penalty_name = "SCAD",
    map = map, n_coef = q,
    params = c("lambda", "a"),
    params_bounds = list(lambda = c(0, Inf), a = c(2, Inf)),
    link_params = list(lambda = link_lambda, a = link_a),
    params_smooth = c(lambda = TRUE, a = TRUE)
  )
}

#' @rdname scad_penalty
#' @export
mcp_penalty <- function(map = NULL, n_coef = 1L,
                        link_lambda = linkfunctions7::log_link(),
                        link_gamma = linkfunctions7::bounded_link(lwr = 1)) {
  q <- if (is.null(map)) as.integer(n_coef) else ncol(map <- as.matrix(map))
  McpPenalty(
    penalty_name = "MCP",
    map = map, n_coef = q,
    params = c("lambda", "gamma"),
    params_bounds = list(lambda = c(0, Inf), gamma = c(1, Inf)),
    link_params = list(lambda = link_lambda, gamma = link_gamma),
    params_smooth = c(lambda = TRUE, gamma = TRUE)
  )
}

#' The Piecewise Regions of SCAD and MCP
#' @description The sign, the absolute value and the region indicators of
#'   each element, shared by every method below.
#' @param t The mapped coefficients \eqn{D\beta}.
#' @param lam,a,gam The hyperparameters.
#' @return A list of vectors.
#' @keywords internal
scad_parts <- function(t, lam, a) {
  s <- sign(t)
  u <- abs(t)
  r1 <- u <= lam
  r3 <- u > a * lam
  r2 <- !r1 & !r3
  list(s = s, u = u, r1 = r1, r2 = r2, r3 = r3)
}

#' @title SCAD Methods
#' @name penalty_value.ScadPenalty
#' @description The closed piecewise forms; see \code{\link{scad_penalty}}.
#' @param pen A \code{ScadPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A list containing \code{lambda} and \code{a}.
#' @param scale Handled by the generic.
#' @param ... Unused.
#' @return See the generic pages.
#' @keywords internal
S7::method(penalty_value, ScadPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; a <- theta$a
  p <- scad_parts(map_apply(pen, beta), lam, a)
  u <- p$u
  sum(ifelse(p$r1, lam * u,
      ifelse(p$r2, (2 * a * lam * u - u^2 - lam^2) / (2 * (a - 1)),
             lam^2 * (a + 1) / 2)))
}

#' @rdname penalty_value.ScadPenalty
#' @name penalty_gradient.ScadPenalty
#' @keywords internal
S7::method(penalty_gradient, ScadPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; a <- theta$a
  p <- scad_parts(map_apply(pen, beta), lam, a)
  d1 <- ifelse(p$r1, lam, ifelse(p$r2, (a * lam - p$u) / (a - 1), 0))
  map_back(pen, p$s * d1)
}

#' @rdname penalty_value.ScadPenalty
#' @name penalty_hessian.ScadPenalty
#' @keywords internal
S7::method(penalty_hessian, ScadPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; a <- theta$a
  p <- scad_parts(map_apply(pen, beta), lam, a)
  map_quad(pen, ifelse(p$r2, -1 / (a - 1), 0))
}

#' @rdname penalty_value.ScadPenalty
#' @name penalty_grad_theta.ScadPenalty
#' @keywords internal
S7::method(penalty_grad_theta, ScadPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    lam <- theta$lambda; a <- theta$a
    p <- scad_parts(map_apply(pen, beta), lam, a)
    u <- p$u
    dl <- ifelse(p$r1, u,
          ifelse(p$r2, (a * u - lam) / (a - 1), lam * (a + 1)))
    da <- ifelse(p$r2, (u - lam)^2 / (2 * (a - 1)^2),
          ifelse(p$r3, lam^2 / 2, 0))
    list(lambda = sum(dl), a = sum(da))
  }

#' @rdname penalty_value.ScadPenalty
#' @name penalty_hess_theta.ScadPenalty
#' @keywords internal
S7::method(penalty_hess_theta, ScadPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    lam <- theta$lambda; a <- theta$a
    p <- scad_parts(map_apply(pen, beta), lam, a)
    u <- p$u
    dll <- ifelse(p$r2, -1 / (a - 1), ifelse(p$r3, a + 1, 0))
    daa <- ifelse(p$r2, -(u - lam)^2 / (a - 1)^3, 0)
    dla <- ifelse(p$r2, (lam - u) / (a - 1)^2, ifelse(p$r3, lam, 0))
    list(lambda_lambda = sum(dll), a_a = sum(daa), lambda_a = sum(dla))
  }

#' @rdname penalty_value.ScadPenalty
#' @name penalty_cross.ScadPenalty
#' @keywords internal
S7::method(penalty_cross, ScadPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    lam <- theta$lambda; a <- theta$a
    p <- scad_parts(map_apply(pen, beta), lam, a)
    dl <- ifelse(p$r1, 1, ifelse(p$r2, a / (a - 1), 0))
    da <- ifelse(p$r2, (p$u - lam) / (a - 1)^2, 0)
    list(lambda = map_back(pen, p$s * dl), a = map_back(pen, p$s * da))
  }

#' @rdname penalty_value.ScadPenalty
#' @name penalty_kinks.ScadPenalty
#' @keywords internal
S7::method(penalty_kinks, ScadPenalty) <- function(pen, theta, ...) {
  lam <- theta$lambda; a <- theta$a
  c(0, -lam, lam, -a * lam, a * lam)
}

#' @rdname penalty_value.ScadPenalty
#' @name is_proper.ScadPenalty
#' @keywords internal
S7::method(is_proper, ScadPenalty) <- function(pen, ...) FALSE

#' @rdname scad_parts
#' @keywords internal
mcp_parts <- function(t, lam, gam) {
  s <- sign(t)
  u <- abs(t)
  r1 <- u <= gam * lam
  list(s = s, u = u, r1 = r1)
}

#' @title MCP Methods
#' @name penalty_value.McpPenalty
#' @description The closed piecewise forms; see \code{\link{mcp_penalty}}.
#' @param pen An \code{McpPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A list containing \code{lambda} and \code{gamma}.
#' @param scale Handled by the generic.
#' @param ... Unused.
#' @return See the generic pages.
#' @keywords internal
S7::method(penalty_value, McpPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; gam <- theta$gamma
  p <- mcp_parts(map_apply(pen, beta), lam, gam)
  u <- p$u
  sum(ifelse(p$r1, lam * u - u^2 / (2 * gam), gam * lam^2 / 2))
}

#' @rdname penalty_value.McpPenalty
#' @name penalty_gradient.McpPenalty
#' @keywords internal
S7::method(penalty_gradient, McpPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; gam <- theta$gamma
  p <- mcp_parts(map_apply(pen, beta), lam, gam)
  map_back(pen, p$s * ifelse(p$r1, lam - p$u / gam, 0))
}

#' @rdname penalty_value.McpPenalty
#' @name penalty_hessian.McpPenalty
#' @keywords internal
S7::method(penalty_hessian, McpPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; gam <- theta$gamma
  p <- mcp_parts(map_apply(pen, beta), lam, gam)
  map_quad(pen, ifelse(p$r1, -1 / gam, 0))
}

#' @rdname penalty_value.McpPenalty
#' @name penalty_grad_theta.McpPenalty
#' @keywords internal
S7::method(penalty_grad_theta, McpPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    lam <- theta$lambda; gam <- theta$gamma
    p <- mcp_parts(map_apply(pen, beta), lam, gam)
    u <- p$u
    list(lambda = sum(ifelse(p$r1, u, gam * lam)),
         gamma = sum(ifelse(p$r1, u^2 / (2 * gam^2), lam^2 / 2)))
  }

#' @rdname penalty_value.McpPenalty
#' @name penalty_hess_theta.McpPenalty
#' @keywords internal
S7::method(penalty_hess_theta, McpPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    lam <- theta$lambda; gam <- theta$gamma
    p <- mcp_parts(map_apply(pen, beta), lam, gam)
    u <- p$u
    list(lambda_lambda = sum(ifelse(p$r1, 0, gam)),
         gamma_gamma = sum(ifelse(p$r1, -u^2 / gam^3, 0)),
         lambda_gamma = sum(ifelse(p$r1, 0, lam)))
  }

#' @rdname penalty_value.McpPenalty
#' @name penalty_cross.McpPenalty
#' @keywords internal
S7::method(penalty_cross, McpPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    lam <- theta$lambda; gam <- theta$gamma
    p <- mcp_parts(map_apply(pen, beta), lam, gam)
    list(lambda = map_back(pen, p$s * ifelse(p$r1, 1, 0)),
         gamma = map_back(pen, p$s * ifelse(p$r1, p$u / gam^2, 0)))
  }

#' @rdname penalty_value.McpPenalty
#' @name penalty_kinks.McpPenalty
#' @keywords internal
S7::method(penalty_kinks, McpPenalty) <- function(pen, theta, ...) {
  lam <- theta$lambda; gam <- theta$gamma
  c(0, -gam * lam, gam * lam)
}

#' @rdname penalty_value.McpPenalty
#' @name is_proper.McpPenalty
#' @keywords internal
S7::method(is_proper, McpPenalty) <- function(pen, ...) FALSE
