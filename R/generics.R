#' @include penalty_class.R
NULL

#' @title Value of a Penalty
#'
#' @description
#' Returns the scalar \eqn{\rho(D\beta; \theta)}, the amount by which a set of
#' coefficients is penalized. The normalizing constant is included whenever
#' the penalty is proper, so the value is exactly the negative log-density of
#' the prior and can be added to a negative log-likelihood without a further
#' term. [is_proper()] reports whether this is the case. An improper penalty
#' keeps whatever constant its branch defines: a rank-deficient quadratic,
#' additive or structured penalty keeps the constant taken over the range of
#' its matrix (a log pseudo-determinant), a separable penalty under a map that
#' is not square, or is singular, returns \eqn{-\sum_i \log f(t_i; \theta)}
#' without a Jacobian term, and SCAD and MCP return the bare \eqn{\rho}.
#'
#' @details
#' # The normalizing constant
#'
#' Dropping the normalizing constant makes no difference to a fit at a fixed
#' hyperparameter, since it does not depend on \eqn{\beta}. It does make the
#' hyperparameter itself unestimable: with the constant dropped, driving
#' \eqn{\lambda} to zero drives the penalty to zero and the joint maximum runs
#' away. With the constant kept, a proper penalty is a density in \eqn{\beta}
#' for every \eqn{\lambda}, so the joint objective has an interior maximum and
#' a marginal criterion has a mode to expand around. The degrees of freedom of
#' a heavy-tailed prior are estimable for the same reason.
#'
#' # The map
#'
#' \eqn{\rho} is evaluated at \eqn{t = D\beta}, not at \eqn{\beta}. With the
#' map `NULL`, which is the default of every constructor that takes a map,
#' \eqn{t = \beta} and no arithmetic is done. A map of \eqn{m} rows gives
#' \eqn{\rho} an argument of length \eqn{m}, so a second-difference map
#' penalizes curvature and leaves the level and the slope unpenalized.
#'
#' @param pen A [penalty()] object of any branch.
#' @param beta A finite numeric vector of length `pen@n_coef`. Coerced with
#'   `as.numeric()` in the generic, so an integer vector or a one-column
#'   matrix is accepted. A vector of another length signals an error.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`. Reordered
#'   and checked against the penalty's open bounds before dispatch. Pass
#'   `list()` for a penalty with no hyperparameters.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return A single number, which is not `NA` when `beta` is finite and
#'   `theta` lies inside the bounds. The value is defined for a finite `beta`
#'   only: at an infinite coefficient a quadratic penalty returns `NaN`.
#'
#' @examples
#' # A ridge over three coefficients, with the Gaussian constant kept.
#' pen <- quadratic_penalty(diag(3))
#' penalty_value(pen, c(1, 0, -1), list(lambda = 2))
#'
#' # Which is exactly the negative log-density of a N(0, 1/lambda) prior.
#' -sum(stats::dnorm(c(1, 0, -1), sd = 1 / sqrt(2), log = TRUE))
#'
#' # The constant is why the value moves with lambda even at beta = 0.
#' sapply(c(0.5, 1, 2, 8),
#'        function(l) penalty_value(pen, c(0, 0, 0), list(lambda = l)))
#'
#' # Under a second-difference map rho is evaluated at the second differences,
#' # so a straight line contributes nothing to the quadratic part and the
#' # value equals the constant, the same value that the zero vector gives.
#' curve <- quadratic_penalty(diag(2), map = diff(diag(4), differences = 2))
#' penalty_value(curve, c(1, 2, 3, 4), list(lambda = 5))
#' penalty_value(curve, c(0, 0, 0, 0), list(lambda = 5))
#' penalty_value(curve, c(1, 2, 4, 8), list(lambda = 5))
#'
#' @seealso [penalty_gradient()] and [penalty_hessian()] for the coefficient
#'   derivatives, [penalty_grad_theta()] for the hyperparameter ones,
#'   [is_proper()] for whether the constant is there, [penalty_kinks()] for
#'   where the value is not differentiable.
#' @export
penalty_value <- S7::new_generic("penalty_value", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })

#' @title Coefficient Derivatives of a Penalty
#'
#' @description
#' `penalty_gradient()` returns \eqn{\partial\rho/\partial\beta} and
#' `penalty_hessian()` returns \eqn{\partial^2\rho/\partial\beta^2}, both in
#' closed form for every shipped branch and both already carried back through
#' the map. These are the two quantities that a penalized fit adds to the score
#' and to the information at each iteration.
#'
#' @details
#' # Derivatives through the map
#'
#' \eqn{\rho} is a function of \eqn{t = D\beta}, so the chain rule gives
#'
#' \deqn{\frac{\partial\rho}{\partial\beta} = D^\top \frac{\partial\rho}{\partial t},
#'   \qquad
#'   \frac{\partial^2\rho}{\partial\beta^2}
#'     = D^\top \frac{\partial^2\rho}{\partial t^2} D,}
#'
#' and the map is applied here, once, so the result can be added directly to a
#' system in \eqn{\beta}. For a separable penalty with a univariate parent the
#' middle matrix is diagonal and the product is formed without building it;
#' with a multivariate parent the middle matrix is block diagonal.
#'
#' # Derivatives at a kink
#'
#' Where \eqn{\rho} is not differentiable the returned gradient is one element
#' of the subdifferential, chosen by the branch, and the Hessian is the value
#' of the branch's own formula at that point. The important case is a lasso at
#' \eqn{\beta_j = 0}, where the subdifferential is the whole interval
#' \eqn{[-\lambda, \lambda]} and a single number cannot represent it. A block
#' with such a coordinate is handled by [penalty_prox()] rather than by a
#' gradient method, and [penalty_kinks()] returns the points where this
#' happens.
#'
#' @param pen A [penalty()] object of any branch.
#' @param beta A finite numeric vector of length `pen@n_coef`, coerced with
#'   `as.numeric()` in the generic. A vector of another length signals an
#'   error.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return `penalty_gradient()` a numeric vector of length `pen@n_coef`;
#'   `penalty_hessian()` a symmetric matrix of that side. It is a base matrix
#'   even where the map is a \pkg{Matrix}, except for a penalty built by
#'   [quadratic_penalty()] with `blocks > 1`, whose Hessian is a `dgCMatrix`
#'   so that the sparse storage is kept.
#'
#' @examples
#' pen <- quadratic_penalty(diag(2))
#' penalty_gradient(pen, c(1, -1), list(lambda = 3))
#' penalty_hessian(pen, c(1, -1), list(lambda = 3))
#'
#' # A quadratic penalty has a constant Hessian, so the gradient is linear.
#' all.equal(penalty_gradient(pen, c(1, -1), list(lambda = 3)),
#'           drop(penalty_hessian(pen, c(1, -1), list(lambda = 3)) %*% c(1, -1)))
#'
#' # The lasso's gradient is lambda times a sign, and at zero it is one
#' # element of the interval [-lambda, lambda] rather than a derivative.
#' penalty_gradient(lasso_penalty(n_coef = 3), c(2, -0.5, 0),
#'                  list(lambda = 1.5))
#'
#' @seealso [penalty_value()] for the quantity differentiated,
#'   [penalty_grad_theta()] for the hyperparameter derivatives,
#'   [penalty_prox()] for the operator to use where the gradient does not
#'   exist, [check_penalty()] to verify both against a numerical route.
#' @export
penalty_gradient <- S7::new_generic("penalty_gradient", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })

#' @rdname penalty_gradient
#' @export
penalty_hessian <- S7::new_generic("penalty_hessian", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })

#' @title Hyperparameter Derivatives of a Penalty
#'
#' @description
#' The three blocks of derivatives in the hyperparameters.
#' `penalty_grad_theta()` returns \eqn{\partial\rho/\partial\theta_k},
#' `penalty_hess_theta()` returns
#' \eqn{\partial^2\rho/\partial\theta_k\partial\theta_l}, and
#' `penalty_cross()` returns the mixed block
#' \eqn{\partial^2\rho/\partial\beta\,\partial\theta_k}, one coefficient vector
#' per hyperparameter. Each is returned on the parameter scale or on the
#' unconstrained scale, as selected by `scale`.
#'
#' @details
#' # Uses of the three blocks
#'
#' The gradient and the Hessian in \eqn{\theta} are the steps of a joint
#' maximization over coefficients and hyperparameters, and the quantities that
#' are differentiated in the Laplace term of a marginal criterion. The mixed
#' block is the off-diagonal corner of that joint system, and it also enters
#' the implicit function theorem: the penalized mode moves with a
#' hyperparameter as \eqn{\partial\hat\beta/\partial\theta_k = -(H + S)^{-1}
#' \partial^2\rho/\partial\beta\,\partial\theta_k}, so the derivative of a
#' criterion computed at the mode requires it.
#'
#' # The two scales
#'
#' Each hyperparameter carries a link \eqn{g} that maps its open interval
#' onto the whole line. With `scale = "link"` the derivatives are taken in
#' \eqn{\eta = g(\theta)}, the coordinates that an unconstrained optimizer
#' moves. The chain rule is applied in the generic body, so a method always
#' returns the parameter scale and a branch written by a user obtains both
#' scales without further code. Writing \eqn{h = g^{-1}},
#'
#' \deqn{\frac{\partial\rho}{\partial\eta_k}
#'   = \frac{\partial\rho}{\partial\theta_k} h_k'(\eta_k), \qquad
#'   \frac{\partial^2\rho}{\partial\eta_k \partial\eta_l}
#'   = \frac{\partial^2\rho}{\partial\theta_k \partial\theta_l}
#'     h_k'(\eta_k) h_l'(\eta_l)
#'   + \delta_{kl} \frac{\partial\rho}{\partial\theta_k} h_k''(\eta_k),}
#'
#' and the mixed block picks up one factor of \eqn{h_k'} and no second-order
#' term, because a reparametrization of the hyperparameters does not change
#' the coefficient direction.
#'
#' # Names of the elements
#'
#' The gradient and the mixed block are keyed by `pen@params`. The Hessian is
#' keyed by pairs, **diagonals first** in `params` order and then the upper
#' off-diagonal pairs, each joined by an underscore: for an elastic net the
#' keys are `lambda_lambda`, `alpha_alpha`, `lambda_alpha`. Only the upper
#' triangle is returned, since the matrix is symmetric.
#'
#' @param pen A [penalty()] object of any branch.
#' @param beta A finite numeric vector of length `pen@n_coef`, coerced with
#'   `as.numeric()` in the generic. A vector of another length signals an
#'   error.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param scale `"parameter"`, the default, gives derivatives in \eqn{\theta}.
#'   `"link"` gives them in the unconstrained values \eqn{\eta = g(\theta)},
#'   which is the scale on which an optimizer works and on which a
#'   hyperparameter's standard error is built.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return All three return a named list.
#'   `penalty_grad_theta()` has one element per hyperparameter, each a single
#'   number, named by `pen@params`.
#'   `penalty_hess_theta()` has \eqn{p(p+1)/2} elements, each a single number,
#'   keyed by pairs as described above.
#'   `penalty_cross()` has one element per hyperparameter, each a numeric
#'   vector of length `pen@n_coef`.
#'
#' @examples
#' pen <- quadratic_penalty(diag(2))
#' penalty_grad_theta(pen, c(1, -1), list(lambda = 3))
#' penalty_hess_theta(pen, c(1, -1), list(lambda = 3))
#' penalty_cross(pen, c(1, -1), list(lambda = 3))
#'
#' # An elastic net has two hyperparameters, so three Hessian keys.
#' enet <- elasticnet_penalty(n_coef = 3)
#' names(penalty_hess_theta(enet, c(1, -0.4, 0.2),
#'                          list(lambda = 2, alpha = 0.4)))
#'
#' # The link scale is the parameter scale times the chain factor. lambda
#' # rides a log link, so h'(eta) is lambda itself.
#' g_par  <- penalty_grad_theta(pen, c(1, -1), list(lambda = 3))
#' g_link <- penalty_grad_theta(pen, c(1, -1), list(lambda = 3),
#'                              scale = "link")
#' all.equal(g_link$lambda, g_par$lambda * 3)
#'
#' @seealso [penalty_value()] for the quantity differentiated,
#'   [penalty_gradient()] for the coefficient derivatives,
#'   [penalty_dhessian()] for the third-order blocks that a marginal criterion
#'   reads, [check_penalty()] to verify all three against numerical routes.
#' @export
penalty_grad_theta <- S7::new_generic("penalty_grad_theta", "pen",
  function(pen, beta, theta, scale = c("parameter", "link"), ...) {
    scale <- match.arg(scale)
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
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
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
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
    if (length(beta) != pen@n_coef) {
      stop(sprintf("'beta' must have length %d.", pen@n_coef), call. = FALSE)
    }
    out <- S7::S7_dispatch()
    if (scale == "link") out <- ptheta_to_link(pen, theta, cross = out)
    out
  })

#' @title Points Where a Penalty Is Not Smooth
#'
#' @description
#' Returns the values of \eqn{t = D\beta} at which some derivative of
#' \eqn{\rho} is discontinuous, so that a numerical reference straddling one of
#' them would measure the break instead of the formula. [check_penalty()]
#' places its grids clear of these points, and a fitting routine can read them
#' to decide whether a block can be passed to a gradient method.
#'
#' @details
#' # Kinks and the other breaks
#'
#' A kink in the strict sense is a point where the value is continuous and the
#' first derivative jumps: \eqn{t = 0} for the lasso, the elastic net, SCAD and
#' MCP, where the subdifferential opens into an interval and a coefficient can
#' be held exactly at zero. The function returns the larger set of points
#' where any derivative breaks, because a method that uses the second
#' derivative needs those points as well. At \eqn{\lambda = 1}, \eqn{a = 3.7}
#' and \eqn{\gamma = 3}, with unit curvature, the points are:
#'
#' | penalty | returned | what breaks there |
#' |---|---|---|
#' | quadratic, ridge, structured, additive | `numeric(0)` | nothing, the value is smooth |
#' | separable with a smooth parent (Gaussian, Student t) | `numeric(0)` | nothing |
#' | lasso, elastic net | `0` | the first derivative |
#' | SCAD | `0, -1, 1, -3.7, 3.7` | the first at 0, the second at \eqn{\pm\lambda} and \eqn{\pm a\lambda} |
#' | MCP | `0, -3, 3` | the first at 0, the second at \eqn{\pm\gamma\lambda} |
#'
#' SCAD's first derivative is continuous at \eqn{\pm\lambda}, where both
#' pieces give \eqn{\lambda}, while its second derivative jumps from \eqn{0} to
#' \eqn{-1/(a-1)}; at \eqn{\pm a\lambda} the second derivative returns to
#' \eqn{0}. MCP's first derivative is continuous at \eqn{\pm\gamma\lambda},
#' and its second derivative jumps there from \eqn{-1/\gamma} to \eqn{0}.
#'
#' # Dependence on the hyperparameters
#'
#' For a separable penalty the points are fixed when the penalty is built,
#' from the parameters of the parent that [distributions7::fixed()] holds at a
#' value: the lasso has its kink at \eqn{0}, and a Laplace parent with its
#' location fixed at \eqn{2} has its kink at \eqn{2}. For SCAD and MCP only
#' \eqn{0} is fixed. Their other points are at \eqn{\pm\lambda/c_j} and
#' \eqn{\pm a\lambda/c_j} for SCAD and at \eqn{\pm\gamma\lambda/c_j} for MCP,
#' where \eqn{c_j} is the curvature of coordinate \eqn{j} (one by default), so
#' the result depends on `theta` and the argument is required.
#'
#' @param pen A [penalty()] object of any branch.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`. Required
#'   even for a smooth penalty, whose result does not depend on it.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return A numeric vector of the points, unsorted and possibly empty, in the
#'   units of \eqn{t = D\beta}. `numeric(0)` for a penalty that is smooth
#'   everywhere.
#'
#' @examples
#' # Smooth penalties return nothing.
#' penalty_kinks(quadratic_penalty(diag(2)), list(lambda = 1))
#'
#' # The lasso breaks at the origin alone.
#' penalty_kinks(lasso_penalty(), list(lambda = 1))
#'
#' # SCAD returns five points: the kink at zero and the two pairs where its
#' # second derivative changes branch.
#' penalty_kinks(scad_penalty(), list(lambda = 1, a = 3.7))
#'
#' # The points other than zero move with the hyperparameters.
#' penalty_kinks(scad_penalty(), list(lambda = 2, a = 3.7))
#'
#' @seealso [penalty_gradient()] for what the first derivative returns at a
#'   kink, [penalty_prox()] for the operator that steps over one,
#'   [has_prox()] for whether the penalty has such an operator,
#'   [check_penalty()] for the caller that uses this to place a grid.
#' @export
penalty_kinks <- S7::new_generic("penalty_kinks", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

#' @title Is a Penalty a Proper Prior?
#'
#' @description
#' `TRUE` when \eqn{\exp(-\rho)} integrates to one over the coefficients, so
#' the value returned by [penalty_value()] is exactly the negative log-density
#' of a prior in \eqn{\beta} and may be added to a negative log-likelihood as
#' it stands. `FALSE` when it does not. In that case the value keeps the
#' constant that its branch defines, if any: a rank-deficient quadratic,
#' additive or structured penalty keeps its pseudo-determinant term, a
#' separable penalty keeps the parent's own constant without a Jacobian term,
#' and SCAD and MCP carry no constant.
#'
#' @details
#' Three situations make a penalty improper, and they arise from different
#' causes.
#'
#' A **rank-deficient quadratic** puts no cost at all on its null directions,
#' so \eqn{\exp(-\rho)} is flat along them and its integral diverges. This is
#' the ordinary case for a smoothing penalty: second differences over
#' \eqn{q} coefficients leave the level and the slope unpenalized and have
#' rank \eqn{q - 2}. The penalty is still usable, because the likelihood
#' supplies the missing curvature, and its normalizing constant is the log
#' **pseudo**-determinant over the range, which [penalty_logpdet()] returns.
#' The same holds for a rank-deficient additive or structured penalty.
#'
#' A **separable penalty** is a density in \eqn{t = D\beta}. It is proper
#' without a map, and under a square nonsingular map, where the value carries
#' the Jacobian term \eqn{-\log\lvert\det D\rvert} and is a density in
#' \eqn{\beta}. Under a map with fewer rows than columns, or a singular square
#' map (a diagonal map with a zero entry, for instance), \eqn{\exp(-\rho)} is
#' flat along the null space of the map, and the predicate is `FALSE`.
#'
#' **SCAD and MCP** are improper for a third reason. Both are defined by their
#' derivative, which is exactly zero beyond \eqn{a\lambda} and
#' \eqn{\gamma\lambda}, so \eqn{\rho} is constant in every direction far from
#' the origin and \eqn{\exp(-\rho)} is not integrable for any constant. They
#' are not the negative log-density of any prior, and the quantities of a
#' marginal criterion are not available for them.
#'
#' The predicate is a property of the penalty as constructed, so it does not
#' depend on `theta` and does not take it.
#'
#' @param pen A [penalty()] object of any branch.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return A single logical.
#'
#' @examples
#' # A full-rank quadratic is a proper Gaussian prior.
#' is_proper(quadratic_penalty(diag(2)))
#'
#' # Second differences leave the level and the slope unpenalized.
#' is_proper(quadratic_penalty(crossprod(diff(diag(4), differences = 2))))
#' penalty_rank(quadratic_penalty(crossprod(diff(diag(4), differences = 2))))
#'
#' # The lasso is a proper Laplace prior; SCAD and MCP are not densities.
#' c(lasso = is_proper(lasso_penalty()),
#'   scad  = is_proper(scad_penalty()),
#'   mcp   = is_proper(mcp_penalty()))
#'
#' # A lasso on second differences leaves a two-dimensional null space.
#' is_proper(lasso_penalty(map = diff(diag(4), differences = 2)))
#'
#' @seealso [penalty_value()] for the value whose constant this predicate
#'   describes, [penalty_logpdet()] for the pseudo-determinant of a deficient
#'   quadratic, [is_quadratic()] and [has_prox()] for the other two
#'   predicates on which a fitting routine selects a method.
#' @export
is_proper <- S7::new_generic("is_proper", "pen")

#' @title Is a Penalty a Quadratic Form?
#'
#' @description
#' `TRUE` for [quadratic_penalty()] (and [ridge_penalty()], which is built
#' on it), [structured_penalty()] and [additive_penalty()], the branches whose
#' value is the quadratic form \eqn{\tfrac{1}{2}\beta^\top S(\theta)\beta}
#' plus a normalizing constant, with a matrix \eqn{S(\theta)} that
#' [penalty_matrix()] returns. `FALSE` for every other branch. The predicate
#' determines how a marginal likelihood criterion treats the penalty: where it
#' is `TRUE` the penalty's contribution to the criterion is a log
#' pseudo-determinant, and [penalty_matrix()], [penalty_rank()],
#' [penalty_null_basis()] and [penalty_logpdet()] supply the pieces. Where it
#' is `FALSE`, [penalty_matrix()], [penalty_rank()] and [penalty_logpdet()]
#' signal an error, and [penalty_null_basis()] answers only for a separable
#' penalty.
#'
#' @details
#' # The quadratic form
#'
#' A quadratic penalty carries a fixed matrix \eqn{P}, a map \eqn{D} and one
#' smoothing parameter \eqn{\lambda}, and is the negative log-density of the
#' Gaussian prior on \eqn{\beta} with precision \eqn{\lambda D^\top P D}:
#'
#' \deqn{\rho(\beta; \lambda)
#'   = \tfrac{\lambda}{2} (D\beta)^\top P (D\beta)
#'   - \tfrac{1}{2}\log^{+}\lvert \lambda D^\top P D \rvert
#'   + \tfrac{r}{2}\log(2\pi),
#'   \qquad \log^{+}\lvert \lambda D^\top P D \rvert
#'     = r \log \lambda + \log^{+}\lvert D^\top P D \rvert,}
#'
#' with \eqn{\log^{+}} the log pseudo-determinant and
#' \eqn{r = \operatorname{rank}(D^\top P D)}. Without a map \eqn{D} is the
#' identity. A structured penalty has the same form with a \pkg{parameters7}
#' matrix parameter \eqn{S(\theta)} in place of \eqn{\lambda D^\top P D}, so
#' its matrix moves with several hyperparameters and its log-determinant is
#' not linear in any one of them.
#'
#' An additive penalty has the matrix \eqn{\sum_k \lambda_k P_k}. A sum of
#' quadratic forms is a quadratic form, so it has all four quantities. Its
#' matrix moves with one hyperparameter per component, where the plain
#' quadratic branch has a single scale \eqn{\lambda}.
#'
#' @param pen A [penalty()] object of any branch.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return A single logical.
#'
#' @examples
#' # The quadratic and structured branches.
#' is_quadratic(quadratic_penalty(diag(2)))
#' is_quadratic(ridge_penalty())
#' is_quadratic(structured_penalty(
#'   parameters7::log_cholesky(2)))
#'
#' # The additive branch, whose matrix moves with one hyperparameter per
#' # component.
#' add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
#' c(quadratic = is_quadratic(add), rank = penalty_rank(add))
#'
#' # A separable penalty, SCAD and MCP are not quadratic.
#' is_quadratic(lasso_penalty())
#'
#' @section Methods:
#' The method on the base class [penalty()] returns `FALSE`. A quadratic
#' branch registers its own method returning `TRUE`, as
#' [is_quadratic.QuadraticPenalty()], [is_quadratic.StructuredPenalty()] and
#' [is_quadratic.AdditivePenalty()] do.
#'
#' @aliases is_quadratic.penalty
#'
#' @seealso [penalty_matrix()] and the three quantities beside it,
#'   [is_proper()] and [has_prox()] for the other two predicates, and
#'   [quadratic_penalty()], [structured_penalty()] and
#'   [additive_penalty()] for the branches.
#' @export
is_quadratic <- S7::new_generic("is_quadratic", "pen")

S7::method(is_quadratic, penalty) <- function(pen, ...) FALSE

#' @title Quantities of a Quadratic Penalty for a Marginal Criterion
#'
#' @description
#' The four quantities that a REML or marginal likelihood criterion requires
#' from a quadratic penalty. `penalty_matrix()` returns the matrix
#' \eqn{S(\theta)} of the quadratic form, `penalty_rank()` its rank,
#' `penalty_null_basis()` an orthonormal basis of the directions that the
#' penalty leaves unpenalized, and `penalty_logpdet()` the log
#' pseudo-determinant with its first two derivatives in the hyperparameters.
#'
#' @details
#' # Role in a marginal criterion
#'
#' A marginal criterion integrates the coefficients out under the prior that
#' the penalty defines, and a Laplace expansion at the penalized mode leaves
#' \eqn{\tfrac{1}{2}\log^{+}\lvert S(\theta)\rvert -
#' \tfrac{1}{2}\log\lvert H + S(\theta)\rvert}. These generics supply the
#' first term. It is a **pseudo**-determinant, taken over the range of
#' \eqn{S} alone, because a smoothing penalty is rank deficient by design and
#' its ordinary determinant is zero.
#'
#' The rank and the null basis are fixed at construction and do not move with
#' `theta`, which is why neither function takes it. A count of the
#' eigenvalues of the assembled \eqn{S(\theta)} above a tolerance would fall
#' as several smoothing parameters spread apart, whereas the null space of a
#' sum of positive semidefinite matrices is the intersection of their null
#' spaces and does not depend on the weights.
#'
#' # The log pseudo-determinant
#'
#' For the plain quadratic branch \eqn{S = \lambda D^\top P D} and
#'
#' \deqn{\log^{+}\lvert S \rvert
#'   = r\log\lambda + \log^{+}\lvert D^\top P D \rvert,}
#'
#' with \eqn{r = \operatorname{rank}(D^\top P D)}, so the gradient is
#' \eqn{r/\lambda} and the second derivative \eqn{-r/\lambda^2}, both exact
#' and both returned. For the structured branch the value and the derivatives
#' are those of the log-determinant of the \pkg{parameters7} matrix parameter
#' ([parameters7::param_logdet()] and its derivatives). All are on the
#' parameter scale, and there is no `scale` argument.
#'
#' # Branches that supply them
#'
#' The quadratic, structured and additive branches, for which [is_quadratic()]
#' is `TRUE`, supply all four. For a separable penalty ([distrib_penalty()]),
#' `penalty_null_basis()` returns an orthonormal basis of the null space of
#' the map, the directions of \eqn{\beta} on which the penalty does not depend,
#' and the other three signal an error. For SCAD, MCP and any other branch
#' without a method of its own, all four signal an error with a message that
#' names [is_quadratic()].
#'
#' @param pen A [penalty()] object for which [is_quadratic()] is `TRUE`, or
#'   a separable penalty for `penalty_null_basis()`.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same. `penalty_matrix()` and `penalty_logpdet()` only.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return `penalty_matrix()` a symmetric `q x q` matrix, where `q` is
#'   `pen@n_coef`.
#'   `penalty_rank()` a single integer, fixed at construction.
#'   `penalty_null_basis()` a `q x (q - r)` matrix with orthonormal columns
#'   spanning the null space, and a `q x 0` matrix when there is none.
#'   The two matrices are base matrices ordinarily, and `dgCMatrix` objects
#'   for a penalty built by [quadratic_penalty()] with `blocks > 1`.
#'   `penalty_logpdet()` a list of three: `value`, a single number; `grad`, a
#'   named list of one number per hyperparameter; and `hess`, a named list of
#'   one number per hyperparameter pair, keyed diagonals first as
#'   [penalty_hess_theta()] is.
#'
#' @examples
#' # Second differences over four coefficients: rank 2, a two-dimensional
#' # null space holding the level and the slope.
#' pen <- quadratic_penalty(crossprod(diff(diag(4), differences = 2)))
#' penalty_rank(pen)
#' N <- penalty_null_basis(pen)
#' dim(N)
#'
#' # The matrix annihilates that basis exactly.
#' max(abs(penalty_matrix(pen, list(lambda = 2)) %*% N))
#'
#' # And a straight line, which the null space spans, costs nothing.
#' penalty_gradient(pen, c(1, 2, 3, 4), list(lambda = 2))
#'
#' # The log pseudo-determinant is r log(lambda) plus a constant, so its
#' # derivatives are r/lambda and -r/lambda^2.
#' lp <- penalty_logpdet(pen, list(lambda = 2))
#' c(lp$grad$lambda, lp$hess$lambda_lambda)
#' c(penalty_rank(pen) / 2, -penalty_rank(pen) / 4)
#'
#' # A separable penalty under a second-difference map has a null basis,
#' # and no matrix.
#' lsep <- lasso_penalty(map = diff(diag(4), differences = 2))
#' dim(penalty_null_basis(lsep))
#' try(penalty_matrix(lsep, list(lambda = 1)))
#'
#' @section Methods:
#' The four methods on the base class [penalty()] signal an error with a
#' message that names [is_quadratic()], because a penalty that is not
#' quadratic has no matrix, rank, null basis or log pseudo-determinant. The
#' quadratic branches register their own methods, documented at
#' [penalty_matrix.QuadraticPenalty()], [penalty_matrix.StructuredPenalty()]
#' and [penalty_matrix.AdditivePenalty()], and a separable penalty registers
#' [penalty_null_basis.DistribPenalty()].
#'
#' @aliases penalty_matrix.penalty
#' @aliases penalty_rank.penalty
#' @aliases penalty_null_basis.penalty
#' @aliases penalty_logpdet.penalty
#'
#' @seealso [is_quadratic()] for the predicate that determines which
#'   branches supply these, [quadratic_penalty()], [structured_penalty()] and
#'   [additive_penalty()] for those branches, [penalty_hess_theta()] for the
#'   names of `penalty_logpdet()`'s `hess`.
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
  stop(paste0("Only a quadratic penalty (see is_quadratic()) or a separable",
              " penalty has a null basis."),
       call. = FALSE)
}
S7::method(penalty_logpdet, penalty) <- function(pen, theta, ...) {
  stop("Only a quadratic penalty has a log pseudo-determinant; see is_quadratic().",
       call. = FALSE)
}


#' @title Interpretable Quantities of a Penalty's Hyperparameters
#'
#' @description
#' For a penalty whose hyperparameters are coordinates of a chart rather than
#' interpretable quantities, returns the quantities to be reported, together
#' with the Jacobian from those coordinates and the scale on which each
#' quantity's interval is built. Returns `NULL` where the hyperparameters are
#' reported as they stand.
#'
#' @details
#' # Multivariate priors
#'
#' A separable penalty whose parent is a multivariate family carries the free
#' values of a matrix parameter as its hyperparameters: the logarithms of the
#' diagonal of a Cholesky factor and the entries below it, which are not
#' interpretable directly. The interpretable quantities are the standard
#' deviations and the correlations of the effects, which
#' [distributions7::mv_derived()] computes and this generic returns. The same
#' distinction is made by [parameters7::param_readable()] for a matrix
#' parameter.
#'
#' The Jacobian gives a delta-method standard error: with \eqn{V} the
#' variance matrix of the hyperparameters, the reported quantities have
#' variance \eqn{JVJ^\top}. The `transform` element names the scale on which
#' each quantity's confidence interval is built and from which it is mapped
#' back, so that a standard deviation stays positive and a correlation stays
#' inside \eqn{(-1, 1)}.
#'
#' # The base method
#'
#' `NULL` means that the hyperparameters are reported as they stand. This is
#' the result for the quadratic, additive, SCAD and MCP branches and for a
#' separable penalty with a univariate parent, whose hyperparameters (a
#' smoothing parameter, a rate, a scale, a shape) are interpretable on their
#' own scale. A structured penalty also returns `NULL`, although its
#' hyperparameters are the free values of a \pkg{parameters7} matrix
#' parameter.
#'
#' @param pen A [penalty()] object of any branch.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return `NULL`, or a list of four as [distributions7::mv_derived()] returns
#'   them: `value`, a named numeric vector of the reported quantities;
#'   `jacobian`, their derivatives with respect to the hyperparameters, one row
#'   per quantity and one column per hyperparameter; `transform`, a character
#'   vector naming the scale each interval is built on (`"log"`, `"atanh"`,
#'   `"identity"`); and `block`, a character vector grouping the quantities for
#'   printing.
#'
#' @examples
#' # A bivariate Gaussian prior on three pairs of effects. Its hyperparameters
#' # are log-Cholesky coordinates.
#' pen <- distrib_penalty(
#'   distributions7::fixed(distributions7::mvgaussian1_distrib(2),
#'                         mu1 = 0, mu2 = 0), n_coef = 6)
#' pen@params
#'
#' # The reported quantities are two standard deviations and a correlation.
#' r <- penalty_readable(pen, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
#'                                 sigma_L2.1 = 0.5))
#' r$value
#' r$transform
#' dim(r$jacobian)
#'
#' # A smoothing parameter is reported as it stands.
#' penalty_readable(quadratic_penalty(diag(2)), list(lambda = 1))
#'
#' @section Methods:
#' The method on the base class [penalty()] returns `NULL`. A separable
#' penalty registers [penalty_readable.DistribPenalty()], which returns the
#' result of [distributions7::mv_derived()] for a multivariate parent and
#' `NULL` for a univariate one.
#'
#' @aliases penalty_readable.penalty
#'
#' @seealso [penalty_value()], [distributions7::mv_derived()] for the
#'   quantities returned for a multivariate prior,
#'   [parameters7::param_readable()] for the same distinction on a matrix
#'   parameter.
#' @export
penalty_readable <- S7::new_generic("penalty_readable", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

S7::method(penalty_readable, penalty) <- function(pen, theta, ...) NULL
