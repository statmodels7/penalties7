#' @include penalty_class.R
NULL

#' @title Value of a Penalty
#'
#' @description
#' Returns the scalar \eqn{\rho(D\beta; \theta)}, the amount a set of
#' coefficients is penalized by. The normalizing constant is included whenever
#' the penalty is proper, so the value is exactly the negative log-density of
#' the prior and can be added to a negative log-likelihood without a further
#' term. An improper penalty returns the bare \eqn{\rho}; ask [is_proper()]
#' which you have.
#'
#' @details
#' # Why the constant is kept
#'
#' Dropping the normalizing constant makes no difference to a fit at a fixed
#' hyperparameter, since it does not depend on \eqn{\beta}. It makes the
#' hyperparameter itself unestimable: with the constant dropped, driving
#' \eqn{\lambda} to zero drives the penalty to zero and the joint maximum runs
#' away. With it kept, a proper penalty is a density in \eqn{\beta} for every
#' \eqn{\lambda}, so the joint objective has an interior maximum and a marginal
#' criterion has something to expand around. Estimating the degrees of freedom
#' of a heavy-tailed prior needs the constant for the same reason.
#'
#' # The map
#'
#' \eqn{\rho} is evaluated at \eqn{t = D\beta}, not at \eqn{\beta}. With the
#' map `NULL`, the default of every constructor, \eqn{t = \beta} and no
#' arithmetic is done. A map of \eqn{m} rows means \eqn{\rho} sees \eqn{m}
#' numbers, so a second-difference map penalizes curvature and leaves the level
#' and the slope alone.
#'
#' @param pen A [penalty()] object of any branch.
#' @param beta A numeric vector of length `pen@n_coef`. Coerced with
#'   `as.numeric()` in the generic, so an integer vector or a one-column
#'   matrix is accepted.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`. Reordered
#'   and checked against the penalty's open bounds before dispatch. Pass
#'   `list()` for a penalty with no hyperparameters.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A single number. Never `NA` for an argument inside the bounds; a
#'   penalty whose value is infinite at some \eqn{\beta} returns `Inf`.
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
#' # A map sends rho the second differences, so only curvature is charged for.
#' # A straight line contributes nothing to the quadratic part and the value
#' # falls to the constant, the same value the zero vector gives.
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
    S7::S7_dispatch()
  })

#' @title Coefficient Derivatives of a Penalty
#'
#' @description
#' `penalty_gradient()` returns \eqn{\partial\rho/\partial\beta} and
#' `penalty_hessian()` returns \eqn{\partial^2\rho/\partial\beta^2}, both in
#' closed form for every shipped branch and both already carried back through
#' the map. These are the two quantities a penalized fit adds to the score and
#' to the information at each iteration.
#'
#' @details
#' # Through the map
#'
#' \eqn{\rho} is a function of \eqn{t = D\beta}, so the chain rule gives
#'
#' \deqn{\frac{\partial\rho}{\partial\beta} = D^\top \frac{\partial\rho}{\partial t},
#'   \qquad
#'   \frac{\partial^2\rho}{\partial\beta^2}
#'     = D^\top \frac{\partial^2\rho}{\partial t^2} D,}
#'
#' and the map is applied here, once, so a caller adds the result straight into
#' a system in \eqn{\beta}. For a separable penalty the middle matrix is
#' diagonal and the product is formed without building it.
#'
#' # At a kink
#'
#' Where \eqn{\rho} is not differentiable the returned gradient is one element
#' of the subdifferential, chosen by the branch, and the Hessian is whatever
#' the branch's own formula gives there. A lasso at \eqn{\beta_j = 0} is the
#' case that matters: the subdifferential is the whole interval
#' \eqn{[-\lambda, \lambda]} and a single number cannot stand for it. Route the
#' block to [penalty_prox()] instead of to a gradient method, and use
#' [penalty_kinks()] to find out where the difficulty is.
#'
#' @param pen A [penalty()] object of any branch.
#' @param beta A numeric vector of length `pen@n_coef`, coerced with
#'   `as.numeric()` in the generic.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return `penalty_gradient()` a numeric vector of length `pen@n_coef`;
#'   `penalty_hessian()` a symmetric matrix of that side. It is a base matrix
#'   even where the map is a \pkg{Matrix}, and a `dgCMatrix` for a penalty
#'   built by [quadratic_penalty()] with `blocks > 1`, which keeps its sparse
#'   storage deliberately.
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
#' The three blocks a consumer needs to move the hyperparameters rather than
#' the coefficients. `penalty_grad_theta()` returns
#' \eqn{\partial\rho/\partial\theta_k}, `penalty_hess_theta()` returns
#' \eqn{\partial^2\rho/\partial\theta_k\partial\theta_l}, and
#' `penalty_cross()` returns the mixed block
#' \eqn{\partial^2\rho/\partial\beta\,\partial\theta_k}, one coefficient vector
#' per hyperparameter. Each answers on the parameter scale or on the
#' unconstrained scale, chosen by `scale`.
#'
#' @details
#' # What each is for
#'
#' The gradient and the Hessian in \eqn{\theta} are what a joint maximization
#' of coefficients and hyperparameters steps with, and what the Laplace term of
#' a marginal criterion differentiates. The mixed block is the off-diagonal
#' corner of that joint system, and it is also what the implicit function
#' theorem needs: the penalized mode moves with a hyperparameter as
#' \eqn{\partial\hat\beta/\partial\theta_k = -(H + S)^{-1}
#' \partial^2\rho/\partial\beta\,\partial\theta_k}, so a criterion computed at
#' the mode cannot be differentiated without it.
#'
#' # The two scales
#'
#' Each hyperparameter carries a link \eqn{g} mapping its open interval onto
#' the whole line. With `scale = "link"` the derivatives are taken in
#' \eqn{\eta = g(\theta)}, the coordinates an unconstrained optimizer moves.
#' The chain rule is applied in the generic body, so a method always returns
#' the parameter scale and a branch written by a user gets both scales for
#' free. Writing \eqn{h = g^{-1}},
#'
#' \deqn{\frac{\partial\rho}{\partial\eta_k}
#'   = \frac{\partial\rho}{\partial\theta_k} h_k'(\eta_k), \qquad
#'   \frac{\partial^2\rho}{\partial\eta_k \partial\eta_l}
#'   = \frac{\partial^2\rho}{\partial\theta_k \partial\theta_l}
#'     h_k'(\eta_k) h_l'(\eta_l)
#'   + \delta_{kl} \frac{\partial\rho}{\partial\theta_k} h_k''(\eta_k),}
#'
#' with the mixed block picking up one factor of \eqn{h_k'} and no second-order
#' term, a reparametrization of the hyperparameters leaving the coefficient
#' direction alone.
#'
#' # The keying
#'
#' The gradient and the mixed block are keyed by `pen@params`. The Hessian is
#' keyed by pairs, **diagonals first** in `params` order and then the upper
#' off-diagonal pairs, each joined by an underscore: for an elastic net the
#' keys are `lambda_lambda`, `alpha_alpha`, `lambda_alpha`. Only the upper
#' triangle is returned, the matrix being symmetric.
#'
#' @param pen A [penalty()] object of any branch.
#' @param beta A numeric vector of length `pen@n_coef`, coerced with
#'   `as.numeric()` in the generic.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param scale `"parameter"`, the default, gives derivatives in \eqn{\theta}.
#'   `"link"` gives them in the unconstrained values \eqn{\eta = g(\theta)},
#'   which is the scale an optimizer works on and the scale a hyperparameter's
#'   standard error is built on.
#' @param ... Passed to methods. No shipped method reads it.
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
#'   [penalty_dhessian()] for the third-order blocks a marginal criterion
#'   reads, [check_penalty()] to verify all three against numerical routes.
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

#' @title Where a Penalty Stops Being Smooth
#'
#' @description
#' Returns the values of \eqn{t = D\beta} at which some derivative of
#' \eqn{\rho} is discontinuous, so that a numerical reference straddling one of
#' them measures the break and not the formula. [check_penalty()] places its
#' grids clear of these points, and a solver consults them to decide whether a
#' block can go to a gradient method at all.
#'
#' @details
#' # The set is wider than the kinks alone
#'
#' A kink proper is a point where the value is continuous and the first
#' derivative jumps: \eqn{t = 0} for the lasso, the elastic net, SCAD and MCP,
#' where the subdifferential opens into an interval and a coefficient can be
#' held exactly at zero. What is returned is the larger set of points where any
#' derivative breaks, because that is what a consumer of the second derivative
#' needs as well. Measured at \eqn{\lambda = 1, a = 3.7, \gamma = 3}:
#'
#' | penalty | returned | what breaks there |
#' |---|---|---|
#' | quadratic, ridge, structured, additive | `numeric(0)` | nothing, the value is smooth |
#' | lasso, elastic net | `0` | the first derivative |
#' | SCAD | `0, -1, 1, -3.7, 3.7` | the first at 0, the second at \eqn{\pm\lambda} and \eqn{\pm a\lambda} |
#' | MCP | `0, -3, 3` | the first at 0, the second at \eqn{\pm\gamma\lambda} |
#'
#' At \eqn{\pm\lambda} SCAD's derivative is continuous, both branches giving
#' \eqn{\lambda}, while its second derivative jumps from \eqn{0} to
#' \eqn{-1/(a-1)}. The same holds for MCP at \eqn{\pm\gamma\lambda}.
#'
#' # The set moves with the hyperparameters
#'
#' Only \eqn{0} is fixed. SCAD's outer points are at \eqn{\pm\lambda} and
#' \eqn{\pm a\lambda} and MCP's at \eqn{\pm\gamma\lambda}, so the answer
#' depends on `theta` and the argument is not optional. This is also what a
#' hyperparameter path is walked over: the size of the kink, rather than the
#' hyperparameter itself, is the quantity that scales the same way across
#' branches.
#'
#' @param pen A [penalty()] object of any branch.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`. Required
#'   even for a smooth penalty, whose answer does not depend on it.
#' @param ... Passed to methods. No shipped method reads it.
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
#' # Which move with the hyperparameters.
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
#' of a prior and may be added to a negative log-likelihood as it stands.
#' `FALSE` when it does not, and then the value is the bare \eqn{\rho} with no
#' constant attached.
#'
#' @details
#' Two things make a penalty improper, and they are different failures.
#'
#' A **rank-deficient quadratic** puts no cost at all on its null directions,
#' so \eqn{\exp(-\rho)} is flat along them and its integral diverges. This is
#' the ordinary case for a smoothing penalty: second differences over
#' \eqn{q} coefficients leave the level and the slope free and have rank
#' \eqn{q - 2}. The penalty is still usable, because the likelihood supplies the
#' missing curvature; what changes is that the normalizing constant is the log
#' **pseudo**-determinant over the range, which [penalty_logpdet()] returns.
#'
#' **SCAD and MCP** are improper for a different reason: both are defined by
#' their derivative, which is exactly zero beyond \eqn{a\lambda} and
#' \eqn{\gamma\lambda}, so \eqn{\rho} is flat in every direction far from the
#' origin and no constant makes \eqn{\exp(-\rho)} integrable. There is no
#' density behind them at all, and no marginal criterion reaches them.
#'
#' The predicate is a statement about the penalty as constructed, so it does
#' not depend on `theta` and takes none.
#'
#' @param pen A [penalty()] object of any branch.
#' @param ... Passed to methods. No shipped method reads it.
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
#' # The lasso is a proper Laplace prior; SCAD and MCP are densities of
#' # nothing.
#' c(lasso = is_proper(lasso_penalty()),
#'   scad  = is_proper(scad_penalty()),
#'   mcp   = is_proper(mcp_penalty()))
#'
#' @seealso [penalty_value()] for the constant this is about,
#'   [penalty_logpdet()] for the pseudo-determinant a deficient quadratic
#'   uses, [is_quadratic()] and [has_prox()] for the other two predicates a
#'   consumer routes on.
#' @export
is_proper <- S7::new_generic("is_proper", "pen")

#' @title Is a Penalty a Quadratic Form With One Scale?
#'
#' @description
#' `TRUE` for [quadratic_penalty()] and [structured_penalty()], the two
#' branches whose value is \eqn{\tfrac{1}{2}(D\beta)^\top S(\theta) (D\beta)}
#' for a matrix that a consumer can ask for. `FALSE` everywhere else. A
#' marginal likelihood criterion routes on this: where it is `TRUE` the
#' penalty's contribution to the criterion is a log determinant, which
#' [penalty_matrix()], [penalty_rank()], [penalty_null_basis()] and
#' [penalty_logpdet()] supply, and where it is `FALSE` those four reject.
#'
#' @details
#' # The quantity behind it
#'
#' A quadratic penalty carries a fixed matrix \eqn{P} and one smoothing
#' parameter \eqn{\lambda}, and is the negative log-density of the Gaussian
#' prior with precision \eqn{\lambda P} on \eqn{D\beta}:
#'
#' \deqn{\rho(\beta; \lambda)
#'   = \tfrac{\lambda}{2} (D\beta)^\top P (D\beta)
#'   - \tfrac{1}{2}\log^{+}\lvert \lambda P \rvert
#'   + \tfrac{r}{2}\log(2\pi),
#'   \qquad \log^{+}\lvert \lambda P \rvert
#'     = r \log \lambda + \log^{+}\lvert P \rvert,}
#'
#' with \eqn{\log^{+}} the log pseudo-determinant and
#' \eqn{r = \operatorname{rank}(P)}. A structured penalty is the same shape
#' with \eqn{S(\theta)} a \pkg{parameters7} matrix parameter in place of
#' \eqn{\lambda P}, so its matrix moves with several hyperparameters and its
#' log-determinant is not linear in any one of them.
#'
#' # The exception to know about
#'
#' [additive_penalty()] answers `FALSE` and yet supplies `penalty_matrix()`,
#' `penalty_rank()` and `penalty_logpdet()`, rejecting only
#' `penalty_null_basis()`. It is quadratic in the coefficients, its matrix
#' being \eqn{\sum_k \lambda_k P_k}, and it is kept out of the predicate. A
#' consumer that routes on `is_quadratic()` alone will not reach the additive
#' branch's marginal quantities, and [check_penalty()] does not test them.
#'
#' @param pen A [penalty()] object of any branch.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A single logical.
#'
#' @examples
#' # The two branches this is TRUE for.
#' is_quadratic(quadratic_penalty(diag(2)))
#' is_quadratic(ridge_penalty())
#' is_quadratic(structured_penalty(
#'   parameters7::log_cholesky(2, role = "precision")))
#'
#' # And the ones it is FALSE for, including the additive branch, which has a
#' # rank all the same.
#' add <- additive_penalty(list(diag(3), diag(c(1, 1, 0))))
#' is_quadratic(add)
#' penalty_rank(add)
#' is_quadratic(lasso_penalty())
#'
#' @seealso [penalty_matrix()] and the three quantities beside it,
#'   [is_proper()] and [has_prox()] for the other two predicates,
#'   [quadratic_penalty()] and [structured_penalty()] for the branches.
#' @export
is_quadratic <- S7::new_generic("is_quadratic", "pen")

S7::method(is_quadratic, penalty) <- function(pen, ...) FALSE

#' @title The Pieces a Marginal Criterion Consumes
#'
#' @description
#' The four quantities a REML or marginal likelihood criterion needs from a
#' quadratic penalty. `penalty_matrix()` returns the matrix \eqn{S(\theta)} of
#' the quadratic form, `penalty_rank()` its rank, `penalty_null_basis()` an
#' orthonormal basis of the directions it does not penalize, and
#' `penalty_logpdet()` the log pseudo-determinant with its first two derivatives
#' in the hyperparameters.
#'
#' @details
#' # Why the criterion needs them
#'
#' A marginal criterion integrates the coefficients out under the prior the
#' penalty describes, and a Laplace expansion at the penalized mode leaves
#' \eqn{\tfrac{1}{2}\log^{+}\lvert S(\theta)\rvert -
#' \tfrac{1}{2}\log\lvert H + S(\theta)\rvert}. The first term is what these
#' generics supply. It is a **pseudo**-determinant, taken over the range of
#' \eqn{S} alone, because a smoothing penalty is rank deficient by design and
#' the ordinary determinant is zero.
#'
#' The rank and the null basis are fixed at construction and do not move with
#' `theta`, which is why neither takes it. That is deliberate and it matters:
#' counting eigenvalues of the assembled \eqn{S(\theta)} above a tolerance
#' gives a count that falls as several smoothing parameters spread apart, while
#' the null space of a sum of positive semidefinite matrices is the
#' intersection of theirs and does not move at all.
#'
#' # The log pseudo-determinant
#'
#' For the plain quadratic branch \eqn{S = \lambda D'PD} and
#'
#' \deqn{\log^{+}\lvert S \rvert = r\log\lambda + \log^{+}\lvert P \rvert,}
#'
#' so the gradient is \eqn{r/\lambda} and the second derivative
#' \eqn{-r/\lambda^2}, both exact and both returned. For the structured branch
#' the derivatives come from the \pkg{parameters7} parametrization's own
#' log-determinant contract. Both are on the parameter scale; there is no
#' `scale` argument here.
#'
#' # What rejects
#'
#' All four reject for a penalty that is not quadratic, with a message naming
#' [is_quadratic()]. The exception is [additive_penalty()], which answers
#' `is_quadratic()` as `FALSE` and supplies the first, second and fourth of
#' these anyway; `penalty_null_basis()` is the one it rejects.
#'
#' @param pen A [penalty()] object. Must be quadratic, or additive for the
#'   three it supplies.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same. `penalty_matrix()` and `penalty_logpdet()` only.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return `penalty_matrix()` a symmetric `q x q` matrix, where `q` is
#'   `pen@n_coef`.
#'   `penalty_rank()` a single integer, fixed at construction.
#'   `penalty_null_basis()` a `q x (q - r)` matrix with orthonormal columns
#'   spanning the null space, and a `q x 0` matrix when the penalty is full
#'   rank.
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
#' @seealso [is_quadratic()] for the predicate that gates these,
#'   [quadratic_penalty()], [structured_penalty()] and [additive_penalty()]
#'   for the branches that supply them, [penalty_hess_theta()] for the keying
#'   of `penalty_logpdet()`'s `hess`.
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
#' Returns the quantities a reader reads, for a penalty whose hyperparameters
#' are coordinates of a chart and not the quantities themselves, together with
#' the Jacobian from those coordinates and the scale each one's interval
#' belongs on. Returns `NULL` where the hyperparameters already are the
#' quantities, which is every branch but one.
#'
#' @details
#' # The case this exists for
#'
#' A penalty whose prior is a multivariate family carries the free values of a
#' matrix parameter as its hyperparameters: the logarithms of the diagonal of a
#' Cholesky factor and the entries below it. Nobody reads those. What the prior
#' is about is the standard deviations and the correlations of the effects it
#' describes, so [distributions7::mv_derived()] declares them and this generic
#' passes them through. It is the same distinction
#' [parameters7::param_readable()] makes for a matrix parameter.
#'
#' The Jacobian is what a delta-method standard error needs: with \eqn{V} the
#' variance matrix of the hyperparameters, the reported quantities have
#' variance \eqn{JVJ^\top}. The `transform` element says which scale each
#' quantity's confidence interval should be built on and mapped back from, so
#' that a standard deviation stays positive and a correlation stays inside
#' \eqn{(-1, 1)}.
#'
#' # The base method
#'
#' `NULL` says that the hyperparameters are the quantities and a consumer
#' should report them as they stand. That is the answer for every other branch:
#' a smoothing parameter, a rate and a shape are each read on their own scale
#' already.
#'
#' @param pen A [penalty()] object of any branch.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param ... Passed to methods. No shipped method reads it.
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
#' # are log-Cholesky coordinates, which nobody interprets.
#' pen <- distrib_penalty(
#'   distributions7::fixed(distributions7::mvgaussian_distrib(2),
#'                         mu1 = 0, mu2 = 0), n_coef = 6)
#' pen@params
#'
#' # What the prior is about is two standard deviations and a correlation.
#' r <- penalty_readable(pen, list(sigma_log_L1 = 0.2, sigma_log_L2 = -0.1,
#'                                 sigma_L2.1 = 0.5))
#' r$value
#' r$transform
#' dim(r$jacobian)
#'
#' # A smoothing parameter is already the quantity it names, so there is
#' # nothing to derive.
#' penalty_readable(quadratic_penalty(diag(2)), list(lambda = 1))
#'
#' @seealso [penalty_value()], [distributions7::mv_derived()] for the
#'   declaration this passes through, [parameters7::param_readable()] for the
#'   same distinction on a matrix parameter.
#' @export
penalty_readable <- S7::new_generic("penalty_readable", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

S7::method(penalty_readable, penalty) <- function(pen, theta, ...) NULL
