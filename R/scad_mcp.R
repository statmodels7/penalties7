#' @include generics.R
NULL

#' @title S7 Classes for the SCAD and MCP Penalties
#'
#' @description
#' `ScadPenalty` is the class that [scad_penalty()] builds, and `McpPenalty` is
#' the class that [mcp_penalty()] builds. The literature defines both families
#' by the derivative \eqn{\rho'}, and their value is the closed piecewise
#' antiderivative anchored at \eqn{\rho(0) = 0}. Each class adds one property
#' to [penalty()], `curv`, the curvature of the loss in each coordinate of
#' \eqn{D\beta}, and every quantity is computed from the two hyperparameters,
#' the map and `curv`.
#'
#' @details
#' Both penalties are improper. The derivative \eqn{\rho'} is exactly zero
#' beyond \eqn{a\lambda} for SCAD and beyond \eqn{\gamma\lambda} for MCP, so
#' \eqn{\rho} is bounded and flat far from the origin and \eqn{\exp(-\rho)}
#' does not integrate over \eqn{\mathbb{R}^q} at any constant. [is_proper()]
#' returns `FALSE` for both, and [penalty_matrix()] and the three quantities
#' beside it signal an error, so the hyperparameters are chosen along a path
#' or by cross-validation.
#'
#' The two classes differ in the shape parameter and its bound: SCAD's is `a`,
#' on \eqn{(2, \infty)}, and MCP's is `gamma`, on \eqn{(1, \infty)}.
#'
#' @section Methods:
#' The two classes register methods for the same generics. For `ScadPenalty`:
#'
#' - value and coefficient derivatives:
#'   [penalty_value()][penalty_value.ScadPenalty],
#'   [penalty_gradient()][penalty_gradient.ScadPenalty],
#'   [penalty_hessian()][penalty_hessian.ScadPenalty];
#' - hyperparameter derivatives:
#'   [penalty_grad_theta()][penalty_grad_theta.ScadPenalty],
#'   [penalty_hess_theta()][penalty_hess_theta.ScadPenalty],
#'   [penalty_cross()][penalty_cross.ScadPenalty];
#' - points of non-smoothness and kind:
#'   [penalty_kinks()][penalty_kinks.ScadPenalty],
#'   [is_proper()][is_proper.ScadPenalty];
#' - proximal operator: [penalty_prox()][penalty_prox.ScadPenalty],
#'   [penalty_prox_spec()][penalty_prox_spec.ScadPenalty].
#'
#' For `McpPenalty`:
#'
#' - value and coefficient derivatives:
#'   [penalty_value()][penalty_value.McpPenalty],
#'   [penalty_gradient()][penalty_gradient.McpPenalty],
#'   [penalty_hessian()][penalty_hessian.McpPenalty];
#' - hyperparameter derivatives:
#'   [penalty_grad_theta()][penalty_grad_theta.McpPenalty],
#'   [penalty_hess_theta()][penalty_hess_theta.McpPenalty],
#'   [penalty_cross()][penalty_cross.McpPenalty];
#' - points of non-smoothness and kind:
#'   [penalty_kinks()][penalty_kinks.McpPenalty],
#'   [is_proper()][is_proper.McpPenalty];
#' - proximal operator: [penalty_prox()][penalty_prox.McpPenalty],
#'   [penalty_prox_spec()][penalty_prox_spec.McpPenalty].
#'
#' Every other generic dispatches to the method of the base class
#' [penalty()]: [is_quadratic()] returns `FALSE`, and [penalty_matrix()] and
#' the three quantities beside it signal an error.
#'
#' @inheritParams penalty
#' @param curv The curvature of the loss in each coordinate of \eqn{D\beta}, a
#'   numeric vector of positive finite values with one entry per coordinate,
#'   or `numeric(0)`, which stores no curvature and gives every coordinate
#'   unit curvature. The constructors [scad_penalty()] and [mcp_penalty()]
#'   build it with [check_curv()].
#'
#' @return An S7 object of class `ScadPenalty` or `McpPenalty`, inheriting from
#'   [penalty()] and carrying its seven properties and `curv`.
#'
#' @seealso [scad_penalty()] and [mcp_penalty()] for the constructors,
#'   [penalty_value.ScadPenalty()] and [penalty_value.McpPenalty()] for the
#'   values, [penalty_prox()] for the proximal operator used in a fit.
#'
#' @examples
#' S7::S7_inherits(scad_penalty(), ScadPenalty)
#' S7::S7_inherits(mcp_penalty(), McpPenalty)
#'
#' # Both add one property to the base class, the curvature per coordinate.
#' setdiff(names(S7::props(scad_penalty())), names(S7::props(lasso_penalty())))
#'
#' @keywords internal
#' @export
ScadPenalty <- S7::new_class(name = "ScadPenalty", parent = penalty,
                              properties = list(curv = S7::class_numeric))

#' @rdname ScadPenalty
#' @keywords internal
#' @export
McpPenalty <- S7::new_class(name = "McpPenalty", parent = penalty,
                             properties = list(curv = S7::class_numeric))

#' @title Construct the SCAD and MCP Penalties
#'
#' @description
#' The two non-convex selection penalties that shrink small coefficients like
#' the lasso and leave large ones unshrunk. Both are defined by their
#' derivative on \eqn{t \ge 0} and extended evenly, and `scad_penalty()` and
#' `mcp_penalty()` build objects whose value is the closed antiderivative
#' anchored at \eqn{\rho(0) = 0}.
#'
#' @details
#' # The two derivatives
#'
#' SCAD, after Fan and Li (2001), keeps the slope \eqn{\lambda} of the lasso
#' up to \eqn{\lambda}, tapers it linearly to zero over the next
#' \eqn{(a-1)\lambda}, and is zero beyond:
#'
#' \deqn{\rho'(t) = \lambda \;\; (t \le \lambda), \qquad
#'   \rho'(t) = \frac{a\lambda - t}{a - 1} \;\; (\lambda < t \le a\lambda),
#'   \qquad \rho'(t) = 0 \;\; (t > a\lambda).}
#'
#' MCP, after Zhang (2010), tapers from the origin:
#'
#' \deqn{\rho'(t) = \left(\lambda - t/\gamma\right)_{+}.}
#'
#' Integrating gives the values, which saturate at \eqn{(a+1)\lambda^2/2} for
#' SCAD and at \eqn{\gamma\lambda^2/2} for MCP. These are the forms for unit
#' curvature; the section on the curvature of the loss gives the form that
#' `curv` produces.
#'
#' # The shape parameter
#'
#' Near zero both penalties behave like a lasso at rate \eqn{\lambda}, so they
#' threshold and produce exact zeros. Far from zero both are flat, so a large
#' coefficient is estimated without shrinkage and without the bias of the
#' lasso. The shape parameter sets how quickly the transition happens:
#' \eqn{a \to \infty} and \eqn{\gamma \to \infty} recover the lasso, smaller
#' values make the transition more abrupt, and MCP approaches hard
#' thresholding as \eqn{\gamma \to 1}. The values used in the literature are
#' \eqn{a = 3.7} and \eqn{\gamma = 3}. Neither constructor supplies a value,
#' which the caller sets or a path sweeps.
#'
#' # Points of non-smoothness
#'
#' Both penalties have a kink at zero, and both change branch again where the
#' taper begins or ends. [penalty_kinks()] returns all of these points: at
#' unit curvature they are `0`, \eqn{\pm\lambda} and \eqn{\pm a\lambda} for
#' SCAD and `0` and \eqn{\pm\gamma\lambda} for MCP, and with `curv` the rate
#' \eqn{\lambda} in them becomes \eqn{\lambda/c_j}. At the points other than
#' the origin the first derivative is continuous and the second jumps. When
#' `beta` is not supplied, [check_penalty()] draws its coefficients away from
#' the points that [penalty_kinks()] returns.
#'
#' # The curvature of the loss
#'
#' Both penalties are defined on the canonical problem
#' \eqn{\tfrac12(z-\theta)^2 + \rho(\theta)}, whose loss has unit curvature
#' (Fan and Lv 2010), and the value \eqn{a = 3.7} of Fan and Li (2001) is
#' suggested for that problem. A likelihood has curvature \eqn{c_j} in
#' coordinate \eqn{j}, and with `curv` the penalty is
#' \eqn{\sum_j c_j\,\rho(t_j;\ \lambda/c_j)}. The slope at zero stays
#' \eqn{\lambda}, the knee moves to \eqn{a\lambda/c_j} (\eqn{\gamma\lambda/c_j}
#' for MCP), and the ratio of the knee to \eqn{\lambda/c_j} is the shape
#' parameter in every coordinate. Coordinate \eqn{j} saturates at
#' \eqn{(a+1)\lambda^2/(2c_j)} for SCAD and \eqn{\gamma\lambda^2/(2c_j)} for
#' MCP, and its second derivative over the taper is \eqn{-c_j/(a-1)} for SCAD
#' and \eqn{-c_j/\gamma} for MCP. The proximal operator with step \eqn{s} is
#' the unscaled operator at step \eqn{s c_j} and rate \eqn{\lambda/c_j}, so it
#' is single-valued while \eqn{s c_j < a - 1} for SCAD and
#' \eqn{s c_j < \gamma} for MCP. At the step \eqn{s = 1/c_j} these conditions
#' are \eqn{a > 2} and \eqn{\gamma > 1} for any data. By default `curv` is
#' `NULL`, no curvature is stored, and every coordinate has \eqn{c_j = 1}.
#'
#' # Improper penalties
#'
#' Both penalties are improper: \eqn{\rho} is bounded, so \eqn{\exp(-\rho)} is
#' not a density at any constant. [is_proper()] and [is_quadratic()] return
#' `FALSE`, and [penalty_matrix()] and the three quantities beside it signal an
#' error. A marginal criterion therefore cannot be used to choose these
#' hyperparameters, which are chosen along a path or by cross-validation.
#'
#' @param map The matrix \eqn{D}, with one column per coefficient, or `NULL`
#'   (the default) for the identity. Given a map, `n_coef` is taken from its
#'   column count and the argument is ignored. A \pkg{Matrix} object keeps
#'   its storage, and a square map whose off-diagonal entries are all zero is
#'   stored as a [Matrix::Diagonal()], so that a diagonal map has a proximal
#'   operator whatever its class.
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
#' @param curv `NULL` (the default) for unit curvature, or positive finite
#'   numbers, one per coordinate of \eqn{D\beta} or a single value repeated.
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
#' Fan, J. and Lv, J. (2010). A selective overview of variable selection in
#' high dimensional feature space. *Statistica Sinica* **20**, 101-148.
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. *Annals of Statistics* **38**, 894-942.
#'
#' @examples
#' # SCAD at the shape used in the literature. A small coefficient is
#' # penalized at the rate of the lasso, a large one by a constant.
#' pen <- scad_penalty(n_coef = 3)
#' th <- list(lambda = 1, a = 3.7)
#' penalty_value(pen, c(0.5, 2, 5), th)
#'
#' # The value saturates at (a + 1) lambda^2 / 2 per coefficient.
#' penalty_value(scad_penalty(n_coef = 1), 100, th)
#' (3.7 + 1) * 1^2 / 2
#'
#' # The gradient at a large coefficient is exactly zero, so that coefficient
#' # is not shrunk.
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
#' # As the shape grows, the value at 2 approaches the lasso's value 2.
#' sapply(c(3.7, 50, 1000),
#'        function(a) penalty_value(scad_penalty(n_coef = 1), 2,
#'                                  list(lambda = 1, a = a)))
#'
#' # With curv = 2 the slope at zero stays lambda and the value saturates at
#' # (a + 1) lambda^2 / (2 * 2).
#' penalty_value(scad_penalty(n_coef = 1, curv = 2), 100, th)
#' (3.7 + 1) * 1^2 / (2 * 2)
#'
#' @seealso [lasso_penalty()] for the convex penalty that both modify,
#'   [penalty_prox()] for the closed proximal operator, defined for steps
#'   below \eqn{a - 1} for SCAD and below \eqn{\gamma} for MCP,
#'   [penalty_kinks()] for the points where they change branch, and
#'   [is_proper()] for the property that excludes a marginal criterion.
#' @export
scad_penalty <- function(map = NULL, n_coef = 1L,
                         link_lambda = linkfunctions7::log_link(),
                         link_a = linkfunctions7::bounded_link(lwr = 2),
                         curv = NULL) {
  if (!is.null(map)) {
    map <- as_map(map)
    # a map whose off-diagonal entries are all zero is a rescaling whatever
    # its class, and is stored as one, as distrib_penalty() does
    if (!methods::is(map, "diagonalMatrix") && nrow(map) == ncol(map) &&
        isTRUE(Matrix::isDiagonal(map))) {
      map <- Matrix::Diagonal(x = as.numeric(Matrix::diag(map)))
    }
  }
  q <- if (is.null(map)) as.integer(n_coef) else ncol(map)
  ScadPenalty(
    penalty_name = "SCAD",
    map = map, n_coef = q,
    params = c("lambda", "a"),
    params_bounds = list(lambda = c(0, Inf), a = c(2, Inf)),
    link_params = list(lambda = link_lambda, a = link_a),
    params_smooth = c(lambda = TRUE, a = TRUE),
    curv = check_curv(curv, if (is.null(map)) q else nrow(map))
  )
}

#' @rdname scad_penalty
#' @export
mcp_penalty <- function(map = NULL, n_coef = 1L,
                        link_lambda = linkfunctions7::log_link(),
                        link_gamma = linkfunctions7::bounded_link(lwr = 1),
                        curv = NULL) {
  if (!is.null(map)) {
    map <- as_map(map)
    # a map whose off-diagonal entries are all zero is a rescaling whatever
    # its class, and is stored as one, as distrib_penalty() does
    if (!methods::is(map, "diagonalMatrix") && nrow(map) == ncol(map) &&
        isTRUE(Matrix::isDiagonal(map))) {
      map <- Matrix::Diagonal(x = as.numeric(Matrix::diag(map)))
    }
  }
  q <- if (is.null(map)) as.integer(n_coef) else ncol(map)
  McpPenalty(
    penalty_name = "MCP",
    map = map, n_coef = q,
    params = c("lambda", "gamma"),
    params_bounds = list(lambda = c(0, Inf), gamma = c(1, Inf)),
    link_params = list(lambda = link_lambda, gamma = link_gamma),
    params_smooth = c(lambda = TRUE, gamma = TRUE),
    curv = check_curv(curv, if (is.null(map)) q else nrow(map))
  )
}

#' The Curvature of Each Coordinate of a SCAD or MCP Penalty
#'
#' @description
#' `check_curv()` validates the curvature by which a SCAD or MCP penalty is
#' scaled and returns it with one value per coordinate of \eqn{D\beta}.
#' `curv_of()` returns the curvature stored in a penalty.
#'
#' @details
#' SCAD and MCP are defined on the canonical problem
#' \eqn{\tfrac12(z - \theta)^2 + \rho(\theta)}, whose loss has unit
#' curvature. Where the loss has curvature \eqn{c_j} in coordinate \eqn{j},
#' the penalty of that coordinate is \eqn{c_j\,\rho(t_j;\ \lambda/c_j)}: the
#' slope at the origin stays \eqn{\lambda}, and the knee moves to
#' \eqn{a\lambda/c_j} (\eqn{\gamma\lambda/c_j} for MCP), so that the ratio of
#' the knee to \eqn{\lambda/c_j} is the shape parameter in every coordinate.
#'
#' `check_curv()` signals an error for an empty vector, a missing, infinite,
#' zero or negative value, and a length other than one or `n`; a single value
#' is repeated `n` times. `curv_of()` returns the stored curvature, or `n`
#' ones where none is stored, which is the penalty as defined without this
#' scaling. A stored curvature is returned whatever `n` is.
#'
#' @param curv `NULL`, or positive finite numbers, one per coordinate of
#'   \eqn{D\beta} or a single value repeated.
#' @param n The number of coordinates of \eqn{D\beta}.
#' @param pen A [ScadPenalty()] or [McpPenalty()].
#'
#' @return `check_curv()` a numeric vector of length `n`, or `numeric(0)`
#'   for `NULL`; `curv_of()` the stored numeric vector, or a vector of `n`
#'   ones.
#'
#' @keywords internal
check_curv <- function(curv, n) {
  if (is.null(curv)) return(numeric(0))
  curv <- as.numeric(curv)
  if (!length(curv) || anyNA(curv) || any(!is.finite(curv)) ||
      any(curv <= 0)) {
    stop("'curv' must be positive and finite.", call. = FALSE)
  }
  if (length(curv) == 1L) curv <- rep(curv, n)
  if (length(curv) != n) {
    stop(sprintf("'curv' must have %d %s, one per coordinate; it has %d.",
                 n, if (n == 1L) "value" else "values", length(curv)),
         call. = FALSE)
  }
  curv
}

#' @rdname check_curv
#' @keywords internal
curv_of <- function(pen, n) {
  if (length(pen@curv)) pen@curv else rep(1, n)
}

#' The Piecewise Regions of SCAD and MCP
#'
#' @description
#' Split the mapped coefficients into the regions over which the formulas of
#' each family are written, and return the sign and the absolute value with
#' them. The value, the coefficient derivatives and the hyperparameter
#' derivatives of the two classes use these helpers, so each formula is one
#' `ifelse()` over indicators and the comparisons are made once.
#'
#' @details
#' SCAD has three regions and MCP has two, and this is the only difference
#' between the two helpers. Every derivative of either family is even in
#' \eqn{t} up to the sign in front, so the formulas are written in
#' `u = abs(t)` and multiplied by `s = sign(t)` where an odd order needs it.
#'
#' At `t = 0` the sign is `0`, so a coefficient exactly at the kink gets a
#' gradient of zero. This is one element of the subdifferential
#' \eqn{[-\lambda, \lambda]}, and the one that a smooth method can use.
#'
#' @param t The mapped coefficients \eqn{D\beta}, a numeric vector.
#' @param lam The rate of each coordinate, \eqn{\lambda/c_j} with \eqn{c_j}
#'   the curvature of the coordinate: a positive numeric vector as long as
#'   `t`, or a single positive number when every curvature is one.
#' @param a SCAD's shape parameter, a single number above 2. `scad_parts()`
#'   only.
#' @param gam MCP's shape parameter, a single number above 1. `mcp_parts()`
#'   only.
#'
#' @return `scad_parts()` a list of five vectors as long as `t`: `s` the sign,
#'   `u` the absolute value, and the logical indicators `r1`
#'   (\eqn{u \le} `lam`), `r2` (the taper, `lam` \eqn{< u \le a} `lam`) and
#'   `r3` (\eqn{u > a} `lam`).
#'   `mcp_parts()` a list of three vectors: `s`, `u`, and `r1` for
#'   \eqn{u \le \gamma} `lam`, whose complement is the flat region.
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
#' constant is added, since SCAD is improper.
#'
#' @details
#' With unit curvature and \eqn{u = \lvert t\rvert} on each coordinate,
#'
#' \deqn{\rho(u) = \lambda u \;\; (u \le \lambda), \qquad
#'   \rho(u) = \frac{2a\lambda u - u^2 - \lambda^2}{2(a-1)}
#'     \;\; (\lambda < u \le a\lambda), \qquad
#'   \rho(u) = \frac{(a+1)\lambda^2}{2} \;\; (u > a\lambda),}
#'
#' the three pieces meeting continuously at \eqn{\lambda} and \eqn{a\lambda}.
#' With curvature \eqn{c_j} (the `curv` property), coordinate \eqn{j}
#' contributes \eqn{c_j\,\rho(u_j;\ \lambda/c_j)}, the same formula with
#' \eqn{\lambda/c_j} in place of \eqn{\lambda} and multiplied by \eqn{c_j}.
#' The value saturates, which leaves a large coefficient unshrunk and makes
#' the penalty improper.
#'
#' @param pen A [ScadPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A list holding `lambda` and `a`, already aligned and
#'   bound-checked by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number, the sum over coordinates. It is bounded above by
#'   \eqn{\sum_j (a+1)\lambda^2/(2c_j)}, the sum running over the coordinates
#'   of \eqn{D\beta}, which is `pen@n_coef * (a + 1) * lambda^2 / 2` under
#'   the identity map at unit curvature.
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
  t <- map_apply(pen, beta)
  cv <- curv_of(pen, length(t))
  lam <- theta$lambda / cv; a <- theta$a
  p <- scad_parts(t, lam, a)
  u <- p$u
  sum(cv * ifelse(p$r1, lam * u,
      ifelse(p$r2, (2 * a * lam * u - u^2 - lam^2) / (2 * (a - 1)),
             lam^2 * (a + 1) / 2)))
}

#' @title Coefficient Derivatives of a SCAD Penalty
#' @name penalty_gradient.ScadPenalty
#'
#' @description
#' `penalty_gradient()` returns SCAD's defining derivative carried back through
#' the map, and `penalty_hessian()` its second derivative carried back through
#' the map. At unit curvature the second derivative is \eqn{-1/(a-1)} in the
#' tapering region and zero elsewhere.
#'
#' @details
#' With unit curvature, on each coordinate, with \eqn{u = \lvert t\rvert} and
#' \eqn{s} its sign,
#'
#' \deqn{\frac{\partial\rho}{\partial t} = s \cdot
#'   \begin{cases} \lambda & u \le \lambda \\
#'     (a\lambda - u)/(a-1) & \lambda < u \le a\lambda \\
#'     0 & u > a\lambda \end{cases}, \qquad
#'   \frac{\partial^2\rho}{\partial t^2} =
#'   \begin{cases} -1/(a-1) & \lambda < u \le a\lambda \\
#'     0 & \text{otherwise.}\end{cases}}
#'
#' With curvature \eqn{c_j} the regions end at \eqn{\lambda/c_j} and
#' \eqn{a\lambda/c_j}, the first derivative is \eqn{s\lambda} up to
#' \eqn{\lambda/c_j} and \eqn{s\,(a\lambda - c_j u)/(a-1)} over the taper,
#' and the second derivative over the taper is \eqn{-c_j/(a-1)}.
#'
#' The Hessian is **negative** where it is not zero, so the penalty is concave
#' over the taper. A penalized objective stays convex only while the curvature
#' of the likelihood in each coordinate exceeds \eqn{1/(a-1)}
#' (\eqn{c_j/(a-1)} with `curv`), and outside that condition a solver that
#' assumes convexity has no guarantee of convergence.
#'
#' At \eqn{t = 0} the sign is zero, so the gradient returned is `0`. This is
#' one element of the subdifferential \eqn{[-\lambda, \lambda]}, and the one
#' with which a smooth method can take a step. Exact zeros are produced by
#' [penalty_prox()], not by the gradient.
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
  t <- map_apply(pen, beta)
  cv <- curv_of(pen, length(t))
  lam <- theta$lambda / cv; a <- theta$a
  p <- scad_parts(t, lam, a)
  d1 <- ifelse(p$r1, lam, ifelse(p$r2, (a * lam - p$u) / (a - 1), 0))
  map_back(pen, p$s * cv * d1)
}

#' @rdname penalty_gradient.ScadPenalty
#' @name penalty_hessian.ScadPenalty
#' @keywords internal
S7::method(penalty_hessian, ScadPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  cv <- curv_of(pen, length(t))
  lam <- theta$lambda / cv; a <- theta$a
  p <- scad_parts(t, lam, a)
  map_quad(pen, cv * ifelse(p$r2, -1 / (a - 1), 0))
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
#' Each region contributes its own expression, and the sums run over the
#' coordinates that fall in it. At unit curvature, in \eqn{u \le \lambda} the
#' value is \eqn{\lambda u}, so \eqn{\partial\rho/\partial\lambda = u} and the
#' shape does not enter; over the taper
#' \eqn{\partial\rho/\partial\lambda = (au - \lambda)/(a-1)} and
#' \eqn{\partial\rho/\partial a = (u - \lambda)^2/(2(a-1)^2)}; beyond
#' \eqn{a\lambda} the value is \eqn{(a+1)\lambda^2/2}, so
#' \eqn{\partial\rho/\partial\lambda = (a+1)\lambda} and
#' \eqn{\partial\rho/\partial a = \lambda^2/2}, and neither depends on the
#' coefficient. With curvature \eqn{c_j} the coordinate contributes
#' \eqn{c_j\,\rho(u_j;\ \lambda/c_j)}, and each derivative follows by the
#' chain rule through \eqn{\lambda/c_j}.
#'
#' A marginal criterion cannot use these derivatives: a Laplace expansion
#' requires the penalty to be twice differentiable at the mode, and a mode
#' with a coefficient held at zero lies on the kink. The derivatives are used
#' in a joint step that estimates the hyperparameters together with the
#' coefficients, and by [check_penalty()], which compares them with
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
#' # The coordinate at 5, past a * lambda, contributes lambda^2 / 2 to the
#' # derivative in a above and zero to the cross derivative in a.
#' penalty_cross(pen, b, th)$a
#'
#' @seealso [penalty_value.ScadPenalty()] for the quantity differentiated,
#'   [penalty_grad_theta()] for the generic and its two scales,
#'   [check_penalty()], which verifies all three.
#' @keywords internal
S7::method(penalty_grad_theta, ScadPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    cv <- curv_of(pen, length(t))
    lam <- theta$lambda / cv; a <- theta$a
    p <- scad_parts(t, lam, a)
    u <- p$u
    dl <- ifelse(p$r1, u,
          ifelse(p$r2, (a * u - lam) / (a - 1), lam * (a + 1)))
    da <- ifelse(p$r2, (u - lam)^2 / (2 * (a - 1)^2),
          ifelse(p$r3, lam^2 / 2, 0))
    list(lambda = sum(dl), a = sum(cv * da))
  }

#' @rdname penalty_grad_theta.ScadPenalty
#' @name penalty_hess_theta.ScadPenalty
#' @keywords internal
S7::method(penalty_hess_theta, ScadPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    cv <- curv_of(pen, length(t))
    lam <- theta$lambda / cv; a <- theta$a
    p <- scad_parts(t, lam, a)
    u <- p$u
    dll <- ifelse(p$r2, -1 / (a - 1), ifelse(p$r3, a + 1, 0))
    daa <- ifelse(p$r2, -(u - lam)^2 / (a - 1)^3, 0)
    dla <- ifelse(p$r2, (lam - u) / (a - 1)^2, ifelse(p$r3, lam, 0))
    list(lambda_lambda = sum(dll / cv), a_a = sum(cv * daa),
         lambda_a = sum(dla))
  }

#' @rdname penalty_grad_theta.ScadPenalty
#' @name penalty_cross.ScadPenalty
#' @keywords internal
S7::method(penalty_cross, ScadPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    cv <- curv_of(pen, length(t))
    lam <- theta$lambda / cv; a <- theta$a
    p <- scad_parts(t, lam, a)
    dl <- ifelse(p$r1, 1, ifelse(p$r2, a / (a - 1), 0))
    da <- ifelse(p$r2, (p$u - lam) / (a - 1)^2, 0)
    list(lambda = map_back(pen, p$s * dl), a = map_back(pen, p$s * cv * da))
  }

#' @title Smoothness and Kind of a SCAD Penalty
#' @name penalty_kinks.ScadPenalty
#'
#' @description
#' `penalty_kinks()` returns the points where SCAD changes branch, and
#' `is_proper()` returns `FALSE`, since SCAD is bounded and so not a density
#' at any constant.
#'
#' @details
#' At unit curvature the points are `0`, \eqn{\pm\lambda} and
#' \eqn{\pm a\lambda}, and they are of two kinds. At the origin the first
#' derivative jumps from \eqn{-\lambda} to \eqn{\lambda}, which is the kink
#' that produces exact zeros. At the other four the first derivative is
#' continuous and the second jumps between \eqn{0} and \eqn{-1/(a-1)}
#' (\eqn{-c_j/(a-1)} with curvature \eqn{c_j}). All of them are returned
#' because a numerical Hessian across any of them is as wrong as a numerical
#' gradient across the origin.
#'
#' With `curv`, the points \eqn{\pm\lambda/c_j} and \eqn{\pm a\lambda/c_j} are
#' returned for each distinct value of \eqn{\lambda/c_j}, so the vector has
#' more than five points when the curvatures differ. Every point except the
#' origin moves with the hyperparameters, so `theta` is required.
#'
#' @param pen A [ScadPenalty()] object.
#' @param theta A list holding `lambda` and `a`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` a numeric vector of distinct points: `0`, then
#'   \eqn{-\lambda/c_j} and \eqn{\lambda/c_j} over the distinct rates, then
#'   \eqn{-a\lambda/c_j} and \eqn{a\lambda/c_j} over the same rates. At unit
#'   curvature these are the five points `0`, \eqn{-\lambda}, \eqn{\lambda},
#'   \eqn{-a\lambda}, \eqn{a\lambda}, in this order.
#'   `is_proper()` the single logical `FALSE`.
#'
#' @examples
#' penalty_kinks(scad_penalty(), list(lambda = 1, a = 3.7))
#' penalty_kinks(scad_penalty(), list(lambda = 2, a = 3.7))
#'
#' # Two curvatures give two rates and nine points.
#' penalty_kinks(scad_penalty(n_coef = 2, curv = c(1, 2)),
#'               list(lambda = 1, a = 3.7))
#' is_proper(scad_penalty())
#'
#' @seealso [penalty_kinks()] and [is_proper()] for the generics,
#'   [penalty_kinks.McpPenalty()] for the sibling family.
#' @keywords internal
S7::method(penalty_kinks, ScadPenalty) <- function(pen, theta, ...) {
  lam <- unique(theta$lambda / curv_of(pen, 1L)); a <- theta$a
  unique(c(0, -lam, lam, -a * lam, a * lam))
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
#' constant is added, since MCP is improper.
#'
#' @details
#' With unit curvature and \eqn{u = \lvert t\rvert} on each coordinate,
#'
#' \deqn{\rho(u) = \lambda u - \frac{u^2}{2\gamma} \;\;
#'     (u \le \gamma\lambda), \qquad
#'   \rho(u) = \frac{\gamma\lambda^2}{2} \;\; (u > \gamma\lambda),}
#'
#' the two pieces meeting continuously at \eqn{\gamma\lambda}. SCAD keeps the
#' slope of the lasso over an interval before it tapers, while MCP tapers from
#' the origin, so its value is a downward parabola up to \eqn{\gamma\lambda}.
#' With curvature \eqn{c_j} (the `curv` property), coordinate \eqn{j}
#' contributes \eqn{c_j\,\rho(u_j;\ \lambda/c_j)}, the same formula with
#' \eqn{\lambda/c_j} in place of \eqn{\lambda} and multiplied by \eqn{c_j}.
#'
#' @param pen An [McpPenalty()] object.
#' @param beta A numeric vector of length `pen@n_coef`, already coerced by the
#'   generic.
#' @param theta A list holding `lambda` and `gamma`, already aligned and
#'   bound-checked by the generic.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A single number, the sum over coordinates. It is bounded above by
#'   \eqn{\sum_j \gamma\lambda^2/(2c_j)}, the sum running over the coordinates
#'   of \eqn{D\beta}, which is `pen@n_coef * gamma * lambda^2 / 2` under the
#'   identity map at unit curvature.
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
  t <- map_apply(pen, beta)
  cv <- curv_of(pen, length(t))
  lam <- theta$lambda / cv; gam <- theta$gamma
  p <- mcp_parts(t, lam, gam)
  u <- p$u
  sum(cv * ifelse(p$r1, lam * u - u^2 / (2 * gam), gam * lam^2 / 2))
}

#' @title Coefficient Derivatives of an MCP Penalty
#' @name penalty_gradient.McpPenalty
#'
#' @description
#' `penalty_gradient()` returns MCP's defining derivative carried back through
#' the map, and `penalty_hessian()` its second derivative carried back through
#' the map. At unit curvature the second derivative is \eqn{-1/\gamma} inside
#' the taper and zero beyond it.
#'
#' @details
#' With unit curvature, on each coordinate, with \eqn{u = \lvert t\rvert} and
#' \eqn{s} its sign,
#'
#' \deqn{\frac{\partial\rho}{\partial t}
#'     = s\,(\lambda - u/\gamma)_{+}, \qquad
#'   \frac{\partial^2\rho}{\partial t^2} =
#'   \begin{cases} -1/\gamma & u \le \gamma\lambda \\
#'     0 & u > \gamma\lambda. \end{cases}}
#'
#' With curvature \eqn{c_j} the taper ends at \eqn{\gamma\lambda/c_j}, the
#' first derivative is \eqn{s\,(\lambda - c_j u/\gamma)_{+}}, and the second
#' derivative inside the taper is \eqn{-c_j/\gamma}.
#'
#' The Hessian is **negative** wherever the penalty is not constant, which is
#' the non-convexity in the name of MCP, and it is constant there. Among the
#' penalties whose derivative is \eqn{\lambda} at the origin and zero beyond
#' \eqn{\gamma\lambda}, MCP has the smallest maximum concavity, the supremum
#' of \eqn{-\rho''} (Zhang 2010). A penalized objective stays convex only while
#' the curvature of the likelihood in each coordinate exceeds \eqn{1/\gamma}
#' (\eqn{c_j/\gamma} with `curv`).
#'
#' At \eqn{t = 0} the sign is zero, so the gradient returned is `0`, one
#' element of the subdifferential \eqn{[-\lambda, \lambda]}. Exact zeros are
#' produced by [penalty_prox()].
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
  t <- map_apply(pen, beta)
  cv <- curv_of(pen, length(t))
  lam <- theta$lambda / cv; gam <- theta$gamma
  p <- mcp_parts(t, lam, gam)
  map_back(pen, p$s * cv * ifelse(p$r1, lam - p$u / gam, 0))
}

#' @rdname penalty_gradient.McpPenalty
#' @name penalty_hessian.McpPenalty
#' @keywords internal
S7::method(penalty_hessian, McpPenalty) <- function(pen, beta, theta, ...) {
  t <- map_apply(pen, beta)
  cv <- curv_of(pen, length(t))
  lam <- theta$lambda / cv; gam <- theta$gamma
  p <- mcp_parts(t, lam, gam)
  map_quad(pen, cv * ifelse(p$r1, -1 / gam, 0))
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
#' At unit curvature, inside the taper the value is
#' \eqn{\lambda u - u^2/(2\gamma)}, so \eqn{\partial\rho/\partial\lambda = u}
#' and \eqn{\partial\rho/\partial\gamma = u^2/(2\gamma^2)}; beyond it the value
#' is \eqn{\gamma\lambda^2/2} and the two derivatives are \eqn{\gamma\lambda}
#' and \eqn{\lambda^2/2}, and neither depends on the coefficient. The second
#' derivative in \eqn{\lambda} is therefore zero inside the taper, where the
#' value is linear in \eqn{\lambda}. With curvature \eqn{c_j} the coordinate
#' contributes \eqn{c_j\,\rho(u_j;\ \lambda/c_j)}, and each derivative follows
#' by the chain rule through \eqn{\lambda/c_j}.
#'
#' As for SCAD, a marginal criterion cannot use these derivatives: the
#' expansion on which such a criterion rests requires a mode at which the
#' penalty is twice differentiable, and a mode with a coefficient at zero lies
#' on the kink.
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
    t <- map_apply(pen, beta)
    cv <- curv_of(pen, length(t))
    lam <- theta$lambda / cv; gam <- theta$gamma
    p <- mcp_parts(t, lam, gam)
    u <- p$u
    list(lambda = sum(ifelse(p$r1, u, gam * lam)),
         gamma = sum(cv * ifelse(p$r1, u^2 / (2 * gam^2), lam^2 / 2)))
  }

#' @rdname penalty_grad_theta.McpPenalty
#' @name penalty_hess_theta.McpPenalty
#' @keywords internal
S7::method(penalty_hess_theta, McpPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    cv <- curv_of(pen, length(t))
    lam <- theta$lambda / cv; gam <- theta$gamma
    p <- mcp_parts(t, lam, gam)
    u <- p$u
    list(lambda_lambda = sum(ifelse(p$r1, 0, gam) / cv),
         gamma_gamma = sum(cv * ifelse(p$r1, -u^2 / gam^3, 0)),
         lambda_gamma = sum(ifelse(p$r1, 0, lam)))
  }

#' @rdname penalty_grad_theta.McpPenalty
#' @name penalty_cross.McpPenalty
#' @keywords internal
S7::method(penalty_cross, McpPenalty) <-
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    t <- map_apply(pen, beta)
    cv <- curv_of(pen, length(t))
    lam <- theta$lambda / cv; gam <- theta$gamma
    p <- mcp_parts(t, lam, gam)
    list(lambda = map_back(pen, p$s * ifelse(p$r1, 1, 0)),
         gamma = map_back(pen, p$s * cv * ifelse(p$r1, p$u / gam^2, 0)))
  }

#' @title Smoothness and Kind of an MCP Penalty
#' @name penalty_kinks.McpPenalty
#'
#' @description
#' `penalty_kinks()` returns the points where MCP changes branch, and
#' `is_proper()` returns `FALSE`, since MCP is bounded and so not a density at
#' any constant.
#'
#' @details
#' At unit curvature the points are `0` and \eqn{\pm\gamma\lambda}, and they
#' differ in the order of the derivative that jumps. At the origin the first
#' derivative jumps from \eqn{-\lambda} to \eqn{\lambda}, which is the kink
#' that produces exact zeros. At \eqn{\pm\gamma\lambda} the first derivative
#' is continuous, both branches giving zero, and the second jumps from
#' \eqn{-1/\gamma} to \eqn{0} (from \eqn{-c_j/\gamma} with curvature
#' \eqn{c_j}). MCP has one pair of outer points where SCAD has two, because its
#' taper starts at the origin instead of at \eqn{\lambda}.
#'
#' With `curv`, the points \eqn{\pm\gamma\lambda/c_j} are returned for each
#' distinct value of \eqn{\lambda/c_j}, so the vector has more than three
#' points when the curvatures differ. The outer points move with the
#' hyperparameters, so `theta` is required.
#'
#' @param pen An [McpPenalty()] object.
#' @param theta A list holding `lambda` and `gamma`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return `penalty_kinks()` a numeric vector of distinct points: `0`, then
#'   \eqn{-\gamma\lambda/c_j} and \eqn{\gamma\lambda/c_j} over the distinct
#'   rates \eqn{\lambda/c_j}. At unit curvature these are the three points
#'   `0`, \eqn{-\gamma\lambda}, \eqn{\gamma\lambda}, in this order.
#'   `is_proper()` the single logical `FALSE`.
#'
#' @examples
#' penalty_kinks(mcp_penalty(), list(lambda = 1, gamma = 3))
#' penalty_kinks(mcp_penalty(), list(lambda = 0.5, gamma = 3))
#'
#' # Two curvatures give two rates and five points.
#' penalty_kinks(mcp_penalty(n_coef = 2, curv = c(1, 2)),
#'               list(lambda = 1, gamma = 3))
#' is_proper(mcp_penalty())
#'
#' @seealso [penalty_kinks()] and [is_proper()] for the generics,
#'   [penalty_kinks.ScadPenalty()] for the sibling family.
#' @keywords internal
S7::method(penalty_kinks, McpPenalty) <- function(pen, theta, ...) {
  lam <- unique(theta$lambda / curv_of(pen, 1L)); gam <- theta$gamma
  unique(c(0, -gam * lam, gam * lam))
}

#' @rdname penalty_kinks.McpPenalty
#' @name is_proper.McpPenalty
#' @keywords internal
S7::method(is_proper, McpPenalty) <- function(pen, ...) FALSE
