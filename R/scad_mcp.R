#' @include generics.R
NULL

#' @title S7 Classes for the Derivative-Defined Penalties
#'
#' @description
#' `ScadPenalty` is the class [scad_penalty()] builds and `McpPenalty` the one
#' [mcp_penalty()] builds. Both are families the literature defines by
#' \eqn{\rho'} rather than by \eqn{\rho}, and whose value follows from the
#' closed piecewise antiderivative anchored at \eqn{\rho(0) = 0}. Neither adds
#' a property to [penalty()]: everything the branch needs is the two
#' hyperparameters and the map, and every quantity is computed from them.
#'
#' @details
#' Both are improper. \eqn{\rho'} is exactly zero beyond \eqn{a\lambda} for
#' SCAD and beyond \eqn{\gamma\lambda} for MCP, so \eqn{\rho} is bounded and
#' flat far from the origin and \eqn{\exp(-\rho)} does not integrate over
#' \eqn{\mathbb{R}^q} at any constant. [is_proper()] returns `FALSE` for both,
#' and the four quantities a marginal criterion reads are unavailable: their
#' hyperparameters have to be chosen along a path or by cross-validation.
#'
#' The two differ in which shape parameter they carry and where it is bounded:
#' SCAD's is `a`, on \eqn{(2, \infty)}, and MCP's is `gamma`, on
#' \eqn{(1, \infty)}.
#'
#' @inheritParams penalty
#'
#' @return An S7 object of class `ScadPenalty` or `McpPenalty`, inheriting from
#'   [penalty()] and carrying its seven properties and no others.
#'
#' @seealso [scad_penalty()] and [mcp_penalty()] for the constructors,
#'   [penalty_value.ScadPenalty()] and [penalty_value.McpPenalty()] for what
#'   each computes, [penalty_prox()] for the operator a fit steps with.
#'
#' @examples
#' S7::S7_inherits(scad_penalty(), ScadPenalty)
#' S7::S7_inherits(mcp_penalty(), McpPenalty)
#'
#' # Neither adds a property to the base class.
#' setdiff(names(S7::props(scad_penalty())),
#'         names(S7::props(quadratic_penalty(diag(1)))))
#'
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
#' The two non-convex selection penalties that shrink small coefficients like
#' the lasso and leave large ones alone. Both are defined by their derivative
#' on \eqn{t \ge 0} and extended evenly, so `scad_penalty()` and
#' `mcp_penalty()` build objects whose value is the closed antiderivative
#' anchored at \eqn{\rho(0) = 0}.
#'
#' @details
#' # The two derivatives
#'
#' SCAD, after Fan and Li (2001), holds the lasso's slope up to \eqn{\lambda},
#' tapers it linearly to zero over the next \eqn{(a-1)\lambda}, and stops:
#'
#' \deqn{\rho'(t) = \lambda \;\; (t \le \lambda), \qquad
#'   \rho'(t) = \frac{a\lambda - t}{a - 1} \;\; (\lambda < t \le a\lambda),
#'   \qquad \rho'(t) = 0 \;\; (t > a\lambda).}
#'
#' MCP, after Zhang (2010), starts tapering at once:
#'
#' \deqn{\rho'(t) = \left(\lambda - t/\gamma\right)_{+}.}
#'
#' Integrating gives the values, which saturate at
#' \eqn{(a+1)\lambda^2/2} and \eqn{\gamma\lambda^2/2}. Both agree with a
#' quadrature of the published derivative to `4e-15` or better at every
#' argument, in each of the three regions and at the boundaries between them.
#'
#' # What the shape parameter buys
#'
#' Near zero both behave like a lasso at rate \eqn{\lambda}, so they threshold
#' and produce exact zeros. Far from zero both are flat, so a large coefficient
#' is estimated without shrinkage, so the lasso's bias is gone. The
#' shape parameter says how quickly the transition happens: \eqn{a \to \infty}
#' and \eqn{\gamma \to \infty} recover the lasso, while small values approach
#' hard thresholding. The literature's defaults are \eqn{a = 3.7} and
#' \eqn{\gamma = 3}, and neither constructor supplies one, the value being
#' the caller's to set or a path's to sweep.
#'
#' # Where they are not smooth
#'
#' Both have a kink at zero, and both change branch again where the taper
#' begins or ends. [penalty_kinks()] returns all of them: `0`, \eqn{\pm\lambda}
#' and \eqn{\pm a\lambda} for SCAD, `0` and \eqn{\pm\gamma\lambda} for MCP. At
#' the outer points the first derivative is continuous and the second jumps.
#' [check_penalty()] asks the object and keeps its grids away from all of them.
#'
#' # What they are not
#'
#' Both are improper: \eqn{\rho} is bounded, so \eqn{\exp(-\rho)} is not a
#' density at any constant. [is_proper()] is `FALSE`, [is_quadratic()] is
#' `FALSE`, and [penalty_matrix()] and its three siblings reject. A marginal
#' criterion cannot reach these hyperparameters; a path or cross-validation
#' can.
#'
#' @param map The matrix \eqn{D}, with one column per coefficient, or `NULL`
#'   (the default) for the identity. Given a map, `n_coef` is taken from its
#'   column count and the argument is ignored.
#' @param n_coef The number of coefficients when `map` is `NULL`. A single
#'   whole number, `1L` by default.
#' @param link_lambda The \pkg{linkfunctions7} link carrying \eqn{\lambda} onto
#'   the whole real line. `linkfunctions7::log_link()` by default,
#'   \eqn{\lambda} being positive.
#' @param link_a The link carrying SCAD's shape parameter, bounded below at 2.
#'   `linkfunctions7::bounded_link(lwr = 2)` by default. `scad_penalty()` only.
#' @param link_gamma The link carrying MCP's shape parameter, bounded below at
#'   1. `linkfunctions7::bounded_link(lwr = 1)` by default. `mcp_penalty()`
#'   only.
#'
#' @return `scad_penalty()` a [ScadPenalty()] object with hyperparameters
#'   `lambda` on \eqn{(0, \infty)} and `a` on \eqn{(2, \infty)}.
#'   `mcp_penalty()` an [McpPenalty()] object with `lambda` on
#'   \eqn{(0, \infty)} and `gamma` on \eqn{(1, \infty)}.
#'
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. *Journal of the American Statistical
#' Association* **96**, 1348-1360.
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. *Annals of Statistics* **38**, 894-942.
#'
#' @examples
#' # SCAD at the literature's shape. Small coefficients pay the lasso's rate,
#' # large ones pay a constant.
#' pen <- scad_penalty(n_coef = 3)
#' th <- list(lambda = 1, a = 3.7)
#' penalty_value(pen, c(0.5, 2, 5), th)
#'
#' # The value saturates at (a + 1) lambda^2 / 2 per coefficient.
#' penalty_value(scad_penalty(n_coef = 1), 100, th)
#' (3.7 + 1) * 1^2 / 2
#'
#' # So the gradient of a large coefficient is exactly zero: no shrinkage.
#' penalty_gradient(pen, c(0.5, 2, 5), th)
#'
#' # MCP saturates sooner, at gamma lambda^2 / 2.
#' mcp <- mcp_penalty(n_coef = 3)
#' penalty_value(mcp, c(0.5, 2, 5), list(lambda = 1, gamma = 3))
#' penalty_gradient(mcp, c(0.5, 2, 5), list(lambda = 1, gamma = 3))
#'
#' # Both are bounded, so neither is a density.
#' c(scad = is_proper(pen), mcp = is_proper(mcp))
#'
#' # Letting the shape run recovers the lasso's linear penalty.
#' sapply(c(3.7, 50, 1000),
#'        function(a) penalty_value(scad_penalty(n_coef = 1), 2,
#'                                  list(lambda = 1, a = a)))
#'
#' @seealso [lasso_penalty()] for the convex penalty these improve on,
#'   [penalty_prox()] for the closed operator each has over its convex
#'   region, [penalty_kinks()] for where they change branch,
#'   [is_proper()] for why no marginal criterion reaches them.
#' @export
scad_penalty <- function(map = NULL, n_coef = 1L,
                         link_lambda = linkfunctions7::log_link(),
                         link_a = linkfunctions7::bounded_link(lwr = 2)) {
  q <- if (is.null(map)) as.integer(n_coef) else ncol(map <- as_map(map))
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
  q <- if (is.null(map)) as.integer(n_coef) else ncol(map <- as_map(map))
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
#'
#' @description
#' Split the mapped coefficients into the regions each family's formulas are
#' written over, and return the sign and the absolute value alongside. Shared
#' by every method of the two branches, so a formula reads as one `ifelse()`
#' over indicators instead of recomputing the comparisons.
#'
#' @details
#' SCAD has three regions and MCP has two, which is the whole difference
#' between the two helpers. Every derivative of either family is even in
#' \eqn{t} up to the sign carried out front, so the formulas are written in
#' `u = abs(t)` and multiplied by `s = sign(t)` where an odd order needs it.
#'
#' At `t = 0` the sign is `0`, so a coefficient sitting exactly at the kink
#' gets a gradient of zero: one element of the subdifferential
#' \eqn{[-\lambda, \lambda]}, chosen for being the one a smooth method can
#' use.
#'
#' @param t The mapped coefficients \eqn{D\beta}, a numeric vector.
#' @param lam The hyperparameter \eqn{\lambda}, a single positive number.
#' @param a SCAD's shape parameter, a single number above 2. `scad_parts()`
#'   only.
#' @param gam MCP's shape parameter, a single number above 1. `mcp_parts()`
#'   only.
#'
#' @return `scad_parts()` a list of five vectors as long as `t`: `s` the sign,
#'   `u` the absolute value, and the logical indicators `r1` (\eqn{u \le
#'   \lambda}), `r2` (the taper) and `r3` (\eqn{u > a\lambda}).
#'   `mcp_parts()` a list of three: `s`, `u`, and `r1` for \eqn{u \le
#'   \gamma\lambda}, the complement needing no name.
#'
#' @seealso [scad_penalty()], [penalty_value.ScadPenalty()],
#'   [penalty_value.McpPenalty()]
#'
#' @keywords internal
scad_parts <- function(t, lam, a) {
  s <- sign(t)
  u <- abs(t)
  r1 <- u <= lam
  r3 <- u > a * lam
  r2 <- !r1 & !r3
  list(s = s, u = u, r1 = r1, r2 = r2, r3 = r3)
}

#' @title Value of a SCAD Penalty
#' @name penalty_value.ScadPenalty
#'
#' @description
#' Returns the closed antiderivative of SCAD's defining \eqn{\rho'}, summed
#' over the coordinates of \eqn{D\beta} and anchored at \eqn{\rho(0) = 0}. No
#' constant is added, SCAD being improper.
#'
#' @details
#' With \eqn{u = \lvert t\rvert} on each coordinate,
#'
#' \deqn{\rho(u) = \lambda u \;\; (u \le \lambda), \qquad
#'   \rho(u) = \frac{2a\lambda u - u^2 - \lambda^2}{2(a-1)}
#'     \;\; (\lambda < u \le a\lambda), \qquad
#'   \rho(u) = \frac{(a+1)\lambda^2}{2} \;\; (u > a\lambda),}
#'
#' the three pieces meeting continuously at \eqn{\lambda} and \eqn{a\lambda}.
#' The value saturates, leaving a large coefficient unshrunk and making the
#' penalty improper.
#'
#' @param pen A [ScadPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A list holding `lambda` and `a`, already aligned and
#'   bound-checked by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number, the sum over coordinates. Bounded above by
#'   `pen@n_coef * (a + 1) * lambda^2 / 2`.
#'
#' @examples
#' pen <- scad_penalty(n_coef = 1)
#' th <- list(lambda = 1, a = 3.7)
#'
#' # One coordinate in each region: linear, tapering, flat.
#' sapply(c(0.5, 2, 100), function(t) penalty_value(pen, t, th))
#'
#' # The flat value is (a + 1) lambda^2 / 2.
#' (3.7 + 1) / 2
#'
#' @seealso [scad_penalty()] for the construction,
#'   [penalty_gradient.ScadPenalty()] for the derivative it integrates,
#'   [penalty_value.McpPenalty()] for the sibling family.
#' @keywords internal
S7::method(penalty_value, ScadPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; a <- theta$a
  p <- scad_parts(map_apply(pen, beta), lam, a)
  u <- p$u
  sum(ifelse(p$r1, lam * u,
      ifelse(p$r2, (2 * a * lam * u - u^2 - lam^2) / (2 * (a - 1)),
             lam^2 * (a + 1) / 2)))
}

#' @title Coefficient Derivatives of a SCAD Penalty
#' @name penalty_gradient.ScadPenalty
#'
#' @description
#' `penalty_gradient()` returns SCAD's defining derivative carried back through
#' the map, and `penalty_hessian()` its second derivative, which is
#' \eqn{-1/(a-1)} in the tapering region and zero elsewhere.
#'
#' @details
#' On each coordinate, with \eqn{u = \lvert t\rvert} and \eqn{s} its sign,
#'
#' \deqn{\frac{\partial\rho}{\partial t} = s \cdot
#'   \begin{cases} \lambda & u \le \lambda \\
#'     (a\lambda - u)/(a-1) & \lambda < u \le a\lambda \\
#'     0 & u > a\lambda \end{cases}, \qquad
#'   \frac{\partial^2\rho}{\partial t^2} =
#'   \begin{cases} -1/(a-1) & \lambda < u \le a\lambda \\
#'     0 & \text{otherwise.}\end{cases}}
#'
#' The Hessian is **negative** where it is not zero, which is the non-convexity
#' the family is named for: the penalty bends the wrong way over the taper. A
#' penalized objective stays convex only while the likelihood's own curvature
#' exceeds \eqn{\lambda/(a-1)}, and a solver that assumes convexity has no
#' guarantee here.
#'
#' At \eqn{t = 0} the sign is zero, so the gradient returned is `0`: one
#' element of the subdifferential \eqn{[-\lambda, \lambda]}, and the one a
#' smooth method can take a step with. Selection happens through
#' [penalty_prox()], not through this.
#'
#' @param pen A [ScadPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A list holding `lambda` and `a`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` a numeric vector of length `pen@n_coef`;
#'   `penalty_hessian()` a symmetric `pen@n_coef` by `pen@n_coef` base matrix,
#'   diagonal under the identity map and negative semidefinite always.
#'
#' @examples
#' pen <- scad_penalty(n_coef = 3)
#' th <- list(lambda = 1, a = 3.7)
#'
#' # Linear, tapering, flat: the gradient falls to zero as u grows.
#' penalty_gradient(pen, c(0.5, 2, 5), th)
#'
#' # The Hessian is -1/(a - 1) on the tapering coordinate alone.
#' diag(penalty_hessian(pen, c(0.5, 2, 5), th))
#' -1 / (3.7 - 1)
#'
#' @seealso [penalty_value.ScadPenalty()] for the quantity differentiated,
#'   [penalty_grad_theta.ScadPenalty()] for the hyperparameter blocks,
#'   [penalty_prox()] for the step that produces exact zeros.
#' @keywords internal
S7::method(penalty_gradient, ScadPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; a <- theta$a
  p <- scad_parts(map_apply(pen, beta), lam, a)
  d1 <- ifelse(p$r1, lam, ifelse(p$r2, (a * lam - p$u) / (a - 1), 0))
  map_back(pen, p$s * d1)
}

#' @rdname penalty_gradient.ScadPenalty
#' @name penalty_hessian.ScadPenalty
#' @keywords internal
S7::method(penalty_hessian, ScadPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; a <- theta$a
  p <- scad_parts(map_apply(pen, beta), lam, a)
  map_quad(pen, ifelse(p$r2, -1 / (a - 1), 0))
}

#' @title Hyperparameter Derivatives of a SCAD Penalty
#' @name penalty_grad_theta.ScadPenalty
#'
#' @description
#' The three blocks in \eqn{(\lambda, a)}, each a closed piecewise form summed
#' over the coordinates. `penalty_grad_theta()` returns the two first
#' derivatives, `penalty_hess_theta()` the three second derivatives, and
#' `penalty_cross()` the mixed block, one coefficient vector per
#' hyperparameter.
#'
#' @details
#' Each region contributes its own expression and the sums run over the
#' coordinates that fall in it. In \eqn{u \le \lambda} the value is
#' \eqn{\lambda u}, so \eqn{\partial\rho/\partial\lambda = u} and the shape
#' does not enter; beyond \eqn{a\lambda} the value is
#' \eqn{(a+1)\lambda^2/2}, so \eqn{\partial\rho/\partial\lambda =
#' (a+1)\lambda} and \eqn{\partial\rho/\partial a = \lambda^2/2}, neither
#' depending on the coefficient at all.
#'
#' These derivatives exist, and no marginal criterion can use them: a Laplace
#' expansion needs the penalty to be twice differentiable at the mode, and a
#' mode with a coefficient held at zero sits on the kink. What they serve is a
#' joint step in which the hyperparameters are estimated alongside the
#' coefficients, and [check_penalty()], which compares them against
#' \pkg{numDeriv}.
#'
#' @param pen A [ScadPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A list holding `lambda` and `a`.
#' @param scale Read by the generic, which applies the chain rule onto the
#'   unconstrained scale after dispatch. These methods always return the
#'   parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` a list of two numbers, `lambda` and `a`.
#'   `penalty_hess_theta()` a list of three numbers, `lambda_lambda`, `a_a` and
#'   `lambda_a`.
#'   `penalty_cross()` a list of two numeric vectors of length `pen@n_coef`,
#'   named `lambda` and `a`.
#'
#' @examples
#' pen <- scad_penalty(n_coef = 3)
#' th <- list(lambda = 1, a = 3.7)
#' b <- c(0.5, 2, 5)
#'
#' penalty_grad_theta(pen, b, th)
#' penalty_hess_theta(pen, b, th)
#'
#' # A coordinate past a * lambda contributes lambda^2 / 2 to the shape
#' # derivative and nothing to the mixed block.
#' penalty_cross(pen, b, th)$a
#'
#' @seealso [penalty_value.ScadPenalty()] for the quantity differentiated,
#'   [penalty_grad_theta()] for the generic and its two scales,
#'   [check_penalty()], which verifies all three.
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

#' @rdname penalty_grad_theta.ScadPenalty
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

#' @rdname penalty_grad_theta.ScadPenalty
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

#' @title Smoothness and Kind of a SCAD Penalty
#' @name penalty_kinks.ScadPenalty
#'
#' @description
#' `penalty_kinks()` returns the five points where SCAD changes branch, and
#' `is_proper()` returns `FALSE`, SCAD being bounded and so not a density at
#' any constant.
#'
#' @details
#' The five points are `0`, \eqn{\pm\lambda} and \eqn{\pm a\lambda}, and they
#' are not alike. At the origin the first derivative jumps from
#' \eqn{-\lambda} to \eqn{\lambda}, which is the kink that produces exact
#' zeros. At the other four the first derivative is continuous, both branches
#' agreeing, and the second jumps between \eqn{0} and \eqn{-1/(a-1)}. All five
#' are returned because a numerical Hessian straddling any of them is as wrong
#' as a numerical gradient straddling the origin.
#'
#' Four of the five move with the hyperparameters, so `theta` is not optional.
#'
#' @param pen A [ScadPenalty()] object.
#' @param theta A list holding `lambda` and `a`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` a numeric vector of five points, in the order
#'   `0`, \eqn{-\lambda}, \eqn{\lambda}, \eqn{-a\lambda}, \eqn{a\lambda}.
#'   `is_proper()` the single logical `FALSE`.
#'
#' @examples
#' penalty_kinks(scad_penalty(), list(lambda = 1, a = 3.7))
#' penalty_kinks(scad_penalty(), list(lambda = 2, a = 3.7))
#' is_proper(scad_penalty())
#'
#' @seealso [penalty_kinks()] and [is_proper()] for the generics,
#'   [penalty_kinks.McpPenalty()] for the sibling family's three points.
#' @keywords internal
S7::method(penalty_kinks, ScadPenalty) <- function(pen, theta, ...) {
  lam <- theta$lambda; a <- theta$a
  c(0, -lam, lam, -a * lam, a * lam)
}

#' @rdname penalty_kinks.ScadPenalty
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

#' @title Value of an MCP Penalty
#' @name penalty_value.McpPenalty
#'
#' @description
#' Returns the closed antiderivative of MCP's defining \eqn{\rho'}, summed over
#' the coordinates of \eqn{D\beta} and anchored at \eqn{\rho(0) = 0}. No
#' constant is added, MCP being improper.
#'
#' @details
#' With \eqn{u = \lvert t\rvert} on each coordinate,
#'
#' \deqn{\rho(u) = \lambda u - \frac{u^2}{2\gamma} \;\;
#'     (u \le \gamma\lambda), \qquad
#'   \rho(u) = \frac{\gamma\lambda^2}{2} \;\; (u > \gamma\lambda),}
#'
#' the two pieces meeting continuously at \eqn{\gamma\lambda}. Where SCAD holds
#' the lasso's slope over an interval before tapering, MCP starts tapering at
#' once, so its value is a downward parabola from the origin.
#'
#' @param pen An [McpPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A list holding `lambda` and `gamma`, already aligned and
#'   bound-checked by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number, the sum over coordinates. Bounded above by
#'   `pen@n_coef * gamma * lambda^2 / 2`.
#'
#' @examples
#' pen <- mcp_penalty(n_coef = 1)
#' th <- list(lambda = 1, gamma = 3)
#'
#' # Inside the taper and beyond it.
#' sapply(c(0.5, 2, 100), function(t) penalty_value(pen, t, th))
#'
#' # The flat value is gamma lambda^2 / 2.
#' 3 / 2
#'
#' @seealso [mcp_penalty()] for the construction,
#'   [penalty_gradient.McpPenalty()] for the derivative it integrates,
#'   [penalty_value.ScadPenalty()] for the sibling family.
#' @keywords internal
S7::method(penalty_value, McpPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; gam <- theta$gamma
  p <- mcp_parts(map_apply(pen, beta), lam, gam)
  u <- p$u
  sum(ifelse(p$r1, lam * u - u^2 / (2 * gam), gam * lam^2 / 2))
}

#' @title Coefficient Derivatives of an MCP Penalty
#' @name penalty_gradient.McpPenalty
#'
#' @description
#' `penalty_gradient()` returns MCP's defining derivative carried back through
#' the map, and `penalty_hessian()` its second derivative, which is
#' \eqn{-1/\gamma} inside the taper and zero beyond it.
#'
#' @details
#' On each coordinate, with \eqn{u = \lvert t\rvert} and \eqn{s} its sign,
#'
#' \deqn{\frac{\partial\rho}{\partial t}
#'     = s\,(\lambda - u/\gamma)_{+}, \qquad
#'   \frac{\partial^2\rho}{\partial t^2} =
#'   \begin{cases} -1/\gamma & u \le \gamma\lambda \\
#'     0 & u > \gamma\lambda. \end{cases}}
#'
#' The Hessian is **negative** wherever the penalty is doing anything, which is
#' the non-convexity MCP is named for, and it is constant there: MCP is the
#' penalty with the least curvature among those that threshold at rate
#' \eqn{\lambda} and reach zero bias by \eqn{\gamma\lambda}. A penalized
#' objective stays convex only while the likelihood's own curvature exceeds
#' \eqn{1/\gamma}.
#'
#' At \eqn{t = 0} the sign is zero, so the gradient returned is `0`: one
#' element of the subdifferential \eqn{[-\lambda, \lambda]}. Selection happens
#' through [penalty_prox()].
#'
#' @param pen An [McpPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A list holding `lambda` and `gamma`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_gradient()` a numeric vector of length `pen@n_coef`;
#'   `penalty_hessian()` a symmetric `pen@n_coef` by `pen@n_coef` base matrix,
#'   diagonal under the identity map and negative semidefinite always.
#'
#' @examples
#' pen <- mcp_penalty(n_coef = 3)
#' th <- list(lambda = 1, gamma = 3)
#'
#' # The gradient falls linearly and is exactly zero past gamma * lambda.
#' penalty_gradient(pen, c(0.5, 2, 5), th)
#'
#' # The curvature is a constant -1/gamma inside the taper.
#' diag(penalty_hessian(pen, c(0.5, 2, 5), th))
#'
#' @seealso [penalty_value.McpPenalty()] for the quantity differentiated,
#'   [penalty_grad_theta.McpPenalty()] for the hyperparameter blocks,
#'   [penalty_gradient.ScadPenalty()] for the sibling family.
#' @keywords internal
S7::method(penalty_gradient, McpPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; gam <- theta$gamma
  p <- mcp_parts(map_apply(pen, beta), lam, gam)
  map_back(pen, p$s * ifelse(p$r1, lam - p$u / gam, 0))
}

#' @rdname penalty_gradient.McpPenalty
#' @name penalty_hessian.McpPenalty
#' @keywords internal
S7::method(penalty_hessian, McpPenalty) <- function(pen, beta, theta, ...) {
  lam <- theta$lambda; gam <- theta$gamma
  p <- mcp_parts(map_apply(pen, beta), lam, gam)
  map_quad(pen, ifelse(p$r1, -1 / gam, 0))
}

#' @title Hyperparameter Derivatives of an MCP Penalty
#' @name penalty_grad_theta.McpPenalty
#'
#' @description
#' The three blocks in \eqn{(\lambda, \gamma)}, each a closed piecewise form
#' summed over the coordinates. `penalty_grad_theta()` returns the two first
#' derivatives, `penalty_hess_theta()` the three second derivatives, and
#' `penalty_cross()` the mixed block, one coefficient vector per
#' hyperparameter.
#'
#' @details
#' Inside the taper the value is \eqn{\lambda u - u^2/(2\gamma)}, so
#' \eqn{\partial\rho/\partial\lambda = u} and
#' \eqn{\partial\rho/\partial\gamma = u^2/(2\gamma^2)}; beyond it the value is
#' \eqn{\gamma\lambda^2/2} and the two derivatives are \eqn{\gamma\lambda} and
#' \eqn{\lambda^2/2}, neither depending on the coefficient. The second
#' derivative in \eqn{\lambda} is therefore zero inside the taper, the value
#' being linear in \eqn{\lambda} there.
#'
#' As for SCAD, these derivatives exist and no marginal criterion can use them:
#' the expansion such a criterion rests on needs a mode at which the penalty is
#' twice differentiable, and a mode with a coefficient at zero is on the kink.
#'
#' @param pen An [McpPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`.
#' @param theta A list holding `lambda` and `gamma`.
#' @param scale Read by the generic, which applies the chain rule onto the
#'   unconstrained scale after dispatch. These methods always return the
#'   parameter scale.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_grad_theta()` a list of two numbers, `lambda` and `gamma`.
#'   `penalty_hess_theta()` a list of three numbers, `lambda_lambda`,
#'   `gamma_gamma` and `lambda_gamma`.
#'   `penalty_cross()` a list of two numeric vectors of length `pen@n_coef`,
#'   named `lambda` and `gamma`.
#'
#' @examples
#' pen <- mcp_penalty(n_coef = 3)
#' th <- list(lambda = 1, gamma = 3)
#' b <- c(0.5, 2, 5)
#'
#' penalty_grad_theta(pen, b, th)
#'
#' # The value is linear in lambda inside the taper, so only the coordinate
#' # beyond gamma * lambda contributes to the second derivative.
#' penalty_hess_theta(pen, b, th)$lambda_lambda
#'
#' @seealso [penalty_value.McpPenalty()] for the quantity differentiated,
#'   [penalty_grad_theta()] for the generic and its two scales,
#'   [check_penalty()], which verifies all three.
#' @keywords internal
S7::method(penalty_grad_theta, McpPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    lam <- theta$lambda; gam <- theta$gamma
    p <- mcp_parts(map_apply(pen, beta), lam, gam)
    u <- p$u
    list(lambda = sum(ifelse(p$r1, u, gam * lam)),
         gamma = sum(ifelse(p$r1, u^2 / (2 * gam^2), lam^2 / 2)))
  }

#' @rdname penalty_grad_theta.McpPenalty
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

#' @rdname penalty_grad_theta.McpPenalty
#' @name penalty_cross.McpPenalty
#' @keywords internal
S7::method(penalty_cross, McpPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    lam <- theta$lambda; gam <- theta$gamma
    p <- mcp_parts(map_apply(pen, beta), lam, gam)
    list(lambda = map_back(pen, p$s * ifelse(p$r1, 1, 0)),
         gamma = map_back(pen, p$s * ifelse(p$r1, p$u / gam^2, 0)))
  }

#' @title Smoothness and Kind of an MCP Penalty
#' @name penalty_kinks.McpPenalty
#'
#' @description
#' `penalty_kinks()` returns the three points where MCP changes branch, and
#' `is_proper()` returns `FALSE`, MCP being bounded and so not a density at any
#' constant.
#'
#' @details
#' The three points are `0` and \eqn{\pm\gamma\lambda}, and they differ in
#' order. At the origin the first derivative jumps from \eqn{-\lambda} to
#' \eqn{\lambda}, the kink that produces exact zeros. At
#' \eqn{\pm\gamma\lambda} the first derivative is continuous, both branches
#' giving zero, and the second jumps from \eqn{-1/\gamma} to \eqn{0}. Both
#' outer points move with the hyperparameters, so `theta` is not optional.
#'
#' MCP has one pair of outer points where SCAD has two, its taper starting at
#' the origin instead of at \eqn{\lambda}.
#'
#' @param pen An [McpPenalty()] object.
#' @param theta A list holding `lambda` and `gamma`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` a numeric vector of three points, in the order
#'   `0`, \eqn{-\gamma\lambda}, \eqn{\gamma\lambda}.
#'   `is_proper()` the single logical `FALSE`.
#'
#' @examples
#' penalty_kinks(mcp_penalty(), list(lambda = 1, gamma = 3))
#' penalty_kinks(mcp_penalty(), list(lambda = 0.5, gamma = 3))
#' is_proper(mcp_penalty())
#'
#' @seealso [penalty_kinks()] and [is_proper()] for the generics,
#'   [penalty_kinks.ScadPenalty()] for the sibling family's five points.
#' @keywords internal
S7::method(penalty_kinks, McpPenalty) <- function(pen, theta, ...) {
  lam <- theta$lambda; gam <- theta$gamma
  c(0, -gam * lam, gam * lam)
}

#' @rdname penalty_kinks.McpPenalty
#' @name is_proper.McpPenalty
#' @keywords internal
S7::method(is_proper, McpPenalty) <- function(pen, ...) FALSE
