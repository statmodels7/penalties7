#' @include generics.R quadratic_penalty.R distrib_penalty.R scad_mcp.R
#' @include structured_penalty.R
NULL

#' @title Proximal Operator of a Penalty
#'
#' @description
#' Returns the point that minimizes the penalty plus a quadratic pull towards
#' \eqn{v},
#'
#' \deqn{\mathrm{prox}_{t\rho}(v) = \arg\min_{\beta}
#'   \left\{ \tfrac{1}{2t}\lVert \beta - v \rVert^{2} +
#'   \rho(\beta;\theta) \right\},}
#'
#' A proximal gradient method applies the operator after each gradient step.
#' With it, a penalty with a kink is minimized without differentiating it at
#' the kink, and a coefficient can be set exactly to zero, where a gradient
#' method only brings it close to zero.
#'
#' @details
#' # The operator of each branch
#'
#' | branch | operator |
#' |---|---|
#' | quadratic, structured | one linear solve, \eqn{(I + tS)^{-1}v}, with \eqn{S} the penalty's Hessian |
#' | Gaussian prior centered at zero | the shrinkage \eqn{v/(1 + t/\sigma^2)} |
#' | Laplace prior | the soft threshold \eqn{\mathrm{sign}(v)(\lvert v\rvert - t\lambda)_{+}} |
#' | elastic net | the soft threshold at \eqn{t\lambda\alpha}, followed by division by \eqn{1 + t\lambda(1-\alpha)} |
#' | SCAD, MCP | the piecewise operators of Fan and Li (2001) and Zhang (2010) |
#' | any other separable parent | the stationary point of smallest objective, coordinate by coordinate |
#'
#' Each closed form is the minimizer of the defining objective, and the page of
#' each method gives its formula.
#'
#' # The stationary condition, where there is no closed form
#'
#' A separable penalty whose parent has a differentiable log-density, and which
#' is not one of the closed-form families, is solved coordinate by coordinate
#' from the stationary condition
#'
#' \deqn{\frac{\beta - v}{t} = \ell^{(y)}(\beta),}
#'
#' whose right-hand side is [distributions7::distrib_grad_y()]. For a density
#' that is not log-concave, such as the Student t of [heavy_penalty()], the
#' condition can have more than one root, and the operator is the root of
#' smallest objective. The segment between \eqn{v} and the mode of the parent
#' is scanned for changes of sign, each change is refined by `uniroot()` to a
#' relative accuracy of about `eps^0.75`, and the candidate of smallest
#' objective is returned. For a unimodal parent the global minimizer lies on
#' that segment, so the value returned is the global minimizer. For a parent
#' with more than one mode, the value returned is a stationary point that need
#' not be the global minimizer; [.prox_root()] gives the details. A parent that
#' declares a kink and is not one of the closed-form families is rejected,
#' because the stationary condition does not hold at a kink and the minimizer
#' of such a penalty can lie there.
#'
#' # Under a map
#'
#' The quadratic branch takes any map, because the objective stays
#' quadratic, and a structured penalty has no map. The separable, SCAD and MCP
#' branches need the identity map or a diagonal one. Under a diagonal map with
#' entries \eqn{d_j}, the operator is the operator without the map evaluated at
#' \eqn{d_j v_j} with the step \eqn{t d_j^2}, and then divided by \eqn{d_j};
#' in this way the standardization of a column enters the operator.
#' A general \eqn{D} mixes the coordinates and turns the problem into the
#' generalized-lasso problem, which requires an algorithm of its own, so such
#' a penalty is rejected with a message naming it. A diagonal map is recognized
#' by its \pkg{Matrix} class, as [map_diagonal()] describes.
#'
#' # Calls that signal an error
#'
#' - a penalty with no operator, such as an [additive_penalty()], with a
#'   message naming [has_prox()];
#' - a separable, SCAD or MCP penalty under a map that is not a \pkg{Matrix}
#'   diagonal object with nonzero entries;
#' - a separable penalty whose parent is read in blocks, because the
#'   coordinates of a block do not separate;
#' - a Laplace or elastic-net parent that is not centered at zero, because the
#'   closed form is written for a parent centered at zero;
#' - a separable parent that declares a kink and has no closed form;
#' - a stationary condition that changes sign neither on the segment between
#'   \eqn{v} and the mode of the parent nor in a widened bracket;
#' - SCAD where \eqn{t c_j \ge a - 1} and MCP where \eqn{t c_j \ge \gamma} in
#'   some coordinate, where \eqn{c_j} is the curvature `curv` of the penalty
#'   and \eqn{t} is multiplied by \eqn{d_j^2} under a diagonal map; the
#'   subproblem is then not convex and the operator is set-valued;
#' - a `step` that is not a single positive finite number, or a `v` whose
#'   length is not `pen@n_coef`.
#'
#' @param pen A [penalty()] object. [has_prox()] tests whether it has an
#'   operator.
#' @param v A numeric vector of length `pen@n_coef`, the point of the quadratic
#'   pull. The vector is coerced with `as.numeric()`.
#' @param step The step length \eqn{t}, a single positive finite number. A
#'   vector is rejected; [penalty_prox_spec()] takes one step per coefficient.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector, holding every name in `pen@params`, each a single number.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A numeric vector of the same length as `v`.
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
#' v <- c(2, 0.3, -1.4)
#'
#' # The lasso: a soft threshold, so the middle coordinate becomes exactly 0.
#' penalty_prox(lasso_penalty(n_coef = 3), v, 1, list(lambda = 1))
#'
#' # A ridge: a shrinkage towards zero, which sets no coordinate to zero.
#' penalty_prox(ridge_penalty(n_coef = 3), v, 1, list(lambda = 1.5))
#' v / (1 + 1 * 1.5)
#'
#' # The elastic net applies the soft threshold and then the shrinkage.
#' penalty_prox(elasticnet_penalty(n_coef = 3), v, 1,
#'              list(lambda = 1, alpha = 0.6))
#'
#' # A heavy-tailed prior has no closed form, so the operator is found from
#' # the stationary condition. With sigma = 1, nu = 4 and step 1 it maps 2 to
#' # 1, as the ridge with lambda = 1 does, and it moves the coordinate 10 much
#' # less than that ridge does.
#' v2 <- c(2, 0.3, 10)
#' penalty_prox(heavy_penalty(n_coef = 3), v2, 1, list(sigma = 1, nu = 4))
#' penalty_prox(ridge_penalty(n_coef = 3), v2, 1, list(lambda = 1))
#'
#' # Each operator minimizes the defining objective. The next call finds the
#' # minimizer for the first coordinate of the lasso.
#' rho <- lasso_penalty(n_coef = 3)
#' obj <- function(b1) {
#'   x <- v; x[1] <- b1
#'   0.5 * sum((x - v)^2) / 1 + penalty_value(rho, x, list(lambda = 1))
#' }
#' stats::optimize(obj, c(-5, 5), tol = 1e-12)$minimum
#'
#' # An additive penalty has no operator, and the call signals an error.
#' try(penalty_prox(additive_penalty(list(diag(3), diag(c(1, 1, 0)))),
#'                  v, 1, list(lambda1 = 1, lambda2 = 1)))
#'
#' @section Methods:
#' The method on the base class [penalty()] signals an error naming the penalty
#' and pointing to [has_prox()]; [AdditivePenalty()] inherits it. The classes
#' with an operator of their own are:
#'
#' - [QuadraticPenalty()]: [`penalty_prox()`][penalty_prox.QuadraticPenalty],
#'   one linear solve;
#' - [StructuredPenalty()]:
#'   [`penalty_prox()`][penalty_prox.StructuredPenalty], the same linear
#'   solve;
#' - [DistribPenalty()]: [`penalty_prox()`][penalty_prox.DistribPenalty],
#'   a closed form for the Gaussian, Laplace and elastic-net parents and the
#'   stationary point of smallest objective for any other parent;
#' - [ScadPenalty()]: [`penalty_prox()`][penalty_prox.ScadPenalty], the
#'   piecewise operator of Fan and Li (2001);
#' - [McpPenalty()]: [`penalty_prox()`][penalty_prox.McpPenalty], the
#'   piecewise operator of Zhang (2010).
#'
#' @aliases penalty_prox.penalty
#'
#' @seealso [has_prox()] to test whether a penalty has an operator,
#'   [penalty_prox_spec()] for the same operator as a table that a compiled
#'   loop can read, [penalty_kinks()] for the points where a gradient method
#'   fails, [penalty_gradient()] for the gradient where the penalty is smooth.
#'
#' @export
penalty_prox <- S7::new_generic("penalty_prox", "pen",
  function(pen, v, step, theta, ...) {
    theta <- align_ptheta(pen, theta)
    v <- as.numeric(v)
    if (length(step) != 1L || is.na(step) || !is.finite(step) || step <= 0) {
      stop("'step' must be a single positive finite number.", call. = FALSE)
    }
    if (length(v) != pen@n_coef) {
      stop(sprintf("'v' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })

S7::method(penalty_prox, penalty) <- function(pen, v, step, theta, ...) {
  stop(sprintf(paste0(
    "'%s' has no proximal operator.\n",
    "  A penalty has one when it is a quadratic or structured penalty, a SCAD\n",
    "  or MCP penalty, or a separable penalty under the identity map or a\n",
    "  diagonal map. has_prox() tests this before a call."),
    pen@penalty_name), call. = FALSE)
}

#' @title Does a Penalty Supply a Proximal Operator?
#'
#' @description
#' Returns `TRUE` when [penalty_prox()] has an operator for this penalty, so
#' that a caller can route a block to a proximal method. Returns `FALSE` when
#' the penalty has none, and the block then goes to a smooth method or to a
#' scheme of its own.
#'
#' @details
#' # What it tests
#'
#' The predicate tests four conditions in turn: whether the class registers a
#' [penalty_prox()] method of its own instead of inheriting the error method of
#' the base class; whether the penalty is quadratic or structured, in which
#' case the operator is a linear solve at any map; whether a separable parent
#' is read in blocks, in which case the coordinates do not separate and the
#' result is `FALSE`; and otherwise whether the map is absent or a \pkg{Matrix}
#' diagonal object with nonzero entries, as [map_diagonal()] reads it.
#'
#' The value for each shipped branch is:
#'
#' | penalty | `has_prox()` |
#' |---|---|
#' | quadratic, at any map | `TRUE` |
#' | structured | `TRUE` |
#' | ridge, lasso, elastic net, heavy-tailed, with no map | `TRUE` |
#' | a separable penalty under a diagonal map with nonzero entries | `TRUE` |
#' | SCAD, MCP, with no map or a diagonal map with nonzero entries | `TRUE` |
#' | additive | `FALSE` |
#' | a separable, SCAD or MCP penalty under a general map | `FALSE` |
#' | a separable penalty read in blocks | `FALSE` |
#'
#' # Conditions on the step and on the parent
#'
#' A `TRUE` value means that the penalty has an operator; a call can still
#' signal an error for some steps and some parents. SCAD and MCP return `TRUE`
#' and reject a step for which \eqn{t c_j} reaches \eqn{a - 1} or
#' \eqn{\gamma} in some coordinate (where \eqn{c_j} is the curvature and
#' \eqn{t} is multiplied by \eqn{d_j^2} under a diagonal map), because the
#' subproblem is then not convex and the operator is set-valued. A Laplace or
#' elastic-net parent that is not centered at zero, and a separable parent that
#' declares a kink and has no closed form, return `TRUE`, and the call signals
#' an error.
#'
#' @param pen A [penalty()] object of any branch. Any other object signals an
#'   error.
#'
#' @return A single logical.
#'
#' @examples
#' c(lasso = has_prox(lasso_penalty()),
#'   scad  = has_prox(scad_penalty()),
#'   ridge = has_prox(ridge_penalty()))
#'
#' # The additive branch has no operator.
#' has_prox(additive_penalty(list(diag(3), diag(c(1, 1, 0)))))
#'
#' # For a separable penalty the value depends on the map: TRUE under a
#' # diagonal map, FALSE under a general one.
#' c(diagonal = has_prox(lasso_penalty(map = Matrix::Diagonal(x = c(1, 2, 3)))),
#'   general  = has_prox(lasso_penalty(map = rbind(c(1, -1, 0), c(0, 1, -1)))))
#'
#' # SCAD has an operator, and a step at a - 1 signals an error.
#' has_prox(scad_penalty(n_coef = 1))
#' try(penalty_prox(scad_penalty(n_coef = 1), 2, 2.7, list(lambda = 1, a = 3.7)))
#'
#' @seealso [penalty_prox()] for the operator, [penalty_prox_spec()] for it as
#'   a table, [is_quadratic()] and [is_proper()] for the other two predicates
#'   that a fitting layer uses to route a block.
#'
#' @export
has_prox <- function(pen) {
  if (!S7::S7_inherits(pen, penalty)) {
    stop("'pen' must be a penalty object.", call. = FALSE)
  }
  m <- S7::method(penalty_prox, S7::S7_class(pen))
  owner <- attr(attr(m, "signature")[[1L]], "name")
  if (identical(owner, "penalty")) return(FALSE)
  if (is_quadratic(pen) || S7::S7_inherits(pen, StructuredPenalty)) {
    return(TRUE)
  }
  # A parent read BLOCKWISE has no operator: the operator acts one coordinate
  # at a time and the coordinates of a block do not separate. Asked here
  # rather than discovered at the call, so a fitting layer routes the block to
  # the scheme that can solve it.
  if (S7::S7_inherits(pen, DistribPenalty) && pen@block > 1L) return(FALSE)
  is.null(pen@map) || !is.null(map_diagonal(pen))
}

#' The Diagonal of a Map That Has One
#'
#' @description
#' Returns the diagonal entries of \eqn{D} when the map is a \pkg{Matrix}
#' diagonal object, and `NULL` when there is no map, when the map has another
#' class, or when a diagonal entry is zero or missing.
#'
#' @details
#' A diagonal map rescales each coordinate on its own, so a separable penalty
#' under a diagonal map stays separable. Standardization is such a map:
#' penalizing a column divided by its own spread \eqn{s_j} is penalizing
#' \eqn{\rho(s_j\beta_j)}, so the scaling leaves the design unchanged and the
#' sparsity of a block is kept. A general \eqn{D} mixes the coordinates and
#' turns the problem into the generalized-lasso problem, which requires a
#' different algorithm.
#'
#' The map is recognized by its class, without inspecting its entries: a
#' \pkg{Matrix} diagonal object stores \eqn{q} numbers, while testing a dense
#' matrix for diagonality would cost \eqn{q^2} operations at every call. A base
#' matrix that is diagonal in value therefore gives `NULL`. [distrib_penalty()],
#' [scad_penalty()] and [mcp_penalty()] convert a square map whose
#' off-diagonal entries are all zero to a [Matrix::Diagonal()] when the penalty
#' is built, so a separable, SCAD or MCP penalty built from such a base matrix
#' has a proximal operator.
#'
#' A zero on the diagonal also gives `NULL`, because the transport divides by
#' the diagonal entries.
#'
#' @param pen A [penalty()] object.
#'
#' @return A numeric vector of length `nrow(pen@map)`, or `NULL`. A unit
#'   diagonal object (`diag = "U"`) gives a vector of ones.
#'
#' @seealso [has_prox()], [penalty_prox()], [spec_diag()]
#'
#' @keywords internal
map_diagonal <- function(pen) {
  m <- pen@map
  if (is.null(m)) return(NULL)
  if (!isS4(m) || !methods::is(m, "diagonalMatrix")) return(NULL)
  d <- if (identical(m@diag, "U")) rep(1, nrow(m)) else as.numeric(m@x)
  if (anyNA(d) || any(d == 0)) return(NULL)
  d
}

# The separable branches split by coordinate under the identity map and under
# a DIAGONAL one, which only rescales each coordinate: with u = d b,
#
#     argmin_b (b - v)^2/(2t) + rho(d b)
#       = argmin_u (u - d v)^2/(2 t d^2) + rho(u),   b = u/d
#
# so the operator is the identity-map one read at the scaled point with the
# step scaled by d^2.
# the same penalty with its map removed, which is what the scaled point is
# handed to: the diagonal has already been applied to the point and the step
#' A Penalty With Its Map Removed
#'
#' @description
#' Returns the same penalty with `map` set to `NULL`, which is the object that
#' the diagonal-map transport passes to the formulas written for the identity
#' map.
#'
#' @details
#' The transport applies the diagonal to the point and to the step before the
#' call and divides the result by the diagonal afterwards, so the formula needs
#' the penalty without its map. The other properties are unchanged; in
#' particular `n_coef` is unaffected, because a diagonal map is square.
#'
#' @param pen A [penalty()] object whose map is diagonal.
#'
#' @return The same object with `map` set to `NULL`.
#'
#' @seealso [.prox_scaling()], [spec_diag()], [map_diagonal()]
#'
#' @keywords internal
.undiag <- function(pen) {
  pen@map <- NULL
  pen
}

#' The Diagonal by Which a Proximal Operator Rescales
#'
#' @description
#' Returns the diagonal of the map when the map is a \pkg{Matrix} diagonal
#' object with nonzero entries, `NULL` when there is no map, and signals an
#' error naming the penalty for any other map. The three outcomes correspond
#' to the three cases of a proximal method: the point is rescaled, the
#' operator is applied as it is, or the call stops.
#'
#' @details
#' The identity behind the rescaling is
#'
#' \deqn{\mathrm{prox}_{t\rho(d\,\cdot)}(v) =
#'   \mathrm{prox}_{t d^2 \rho}(d v)/d,}
#'
#' so a caller with a diagonal map multiplies the point by \eqn{d} and the
#' step by \eqn{d^2}, calls the formula written for the identity map, and
#' divides the result by \eqn{d}.
#'
#' The error is signaled here because the formulas written for the identity
#' map would otherwise ignore the map. [has_prox()] returns `FALSE` for the
#' same penalty, so a caller that tests it first does not reach the error.
#'
#' @param pen A [penalty()] object.
#'
#' @return A numeric vector of diagonal entries, or `NULL` when the map is
#'   absent.
#'
#' @section Errors:
#' The message names the penalty and states that its map is not a \pkg{Matrix}
#' diagonal object with nonzero entries, and that the operator is available
#' only under the identity map and such a diagonal map, because under a
#' general map it does not split by coordinate.
#'
#' @seealso [map_diagonal()], [.undiag()], [penalty_prox()]
#'
#' @keywords internal
.prox_scaling <- function(pen) {
  if (is.null(pen@map)) return(NULL)
  d <- map_diagonal(pen)
  if (is.null(d)) {
    stop(sprintf(paste0(
      "'%s' carries a map that is not a Matrix diagonal object with nonzero\n",
      "  entries. The proximal operator is available only under the identity\n",
      "  map and such a diagonal map, because under a general map it\n",
      "  does not split by coordinate."),
      pen@penalty_name), call. = FALSE)
  }
  d
}

#' @title Proximal Operator of a Quadratic or Structured Penalty
#' @name penalty_prox.QuadraticPenalty
#'
#' @description
#' Computes the proximal point with one linear solve, \eqn{(I + tS)^{-1}v},
#' with \eqn{S} the penalty's Hessian: \eqn{\lambda D'PD} for a quadratic
#' penalty and \eqn{\Omega(\theta)} for a structured one. The two methods have
#' the same body, because the objective is quadratic in both.
#'
#' @details
#' The subproblem is
#' \eqn{\arg\min_\beta \{\tfrac{1}{2t}\lVert\beta - v\rVert^2 +
#' \tfrac{1}{2}\beta'S\beta\}}, whose stationary condition is
#' \eqn{(\beta - v)/t + S\beta = 0}. The normalizing constant does not depend
#' on \eqn{\beta} and drops out.
#'
#' The objective stays quadratic under any linear map, so this route places no
#' restriction on \eqn{D}; the separable branches need the identity map or a
#' diagonal one, and a structured penalty has no map. When \eqn{S} is a base
#' matrix the solve is dense and costs \eqn{O(q^3)}. For a penalty built with
#' `blocks > 1`, \eqn{S} is a `dgCMatrix`, the sum
#' `diag(length(v)) + step * S` stays sparse, and `solve()` uses the sparse
#' factorization of \pkg{Matrix}.
#'
#' @param pen A [QuadraticPenalty()] or [StructuredPenalty()] object.
#' @param v A numeric vector of length `pen@n_coef`.
#' @param step The step length \eqn{t}, a single positive finite number.
#' @param theta A named list of hyperparameter values.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the same length as `v`.
#'
#' @examples
#' v <- c(2, 0.3, -1.4)
#'
#' # With S = lambda I the solve is a plain shrinkage.
#' penalty_prox(quadratic_penalty(diag(3)), v, 1, list(lambda = 1.5))
#' v / (1 + 1.5)
#'
#' # A first-difference penalty leaves constant vectors unpenalized, so it
#' # pulls the coefficients towards their mean and not towards zero.
#' penalty_prox(quadratic_penalty(crossprod(diff(diag(3)))), v, 1,
#'              list(lambda = 5))
#'
#' # A structured penalty at a zero log-Cholesky free vector is the ridge with
#' # lambda = 1.
#' s <- structured_penalty(parameters7::log_cholesky(3))
#' penalty_prox(s, v, 1,
#'              as.list(stats::setNames(rep(0, 6), s@params)))
#'
#' @seealso [penalty_prox()] for the generic and the other branches,
#'   [penalty_hessian()] for the matrix against which the system is solved.
#' @keywords internal
S7::method(penalty_prox, QuadraticPenalty) <- function(pen, v, step, theta, ...) {
  S <- penalty_hessian(pen, v, theta)
  as.numeric(solve(diag(length(v)) + step * S, v))
}

#' @rdname penalty_prox.QuadraticPenalty
#' @name penalty_prox.StructuredPenalty
#' @keywords internal
S7::method(penalty_prox, StructuredPenalty) <- function(pen, v, step, theta, ...) {
  S <- penalty_hessian(pen, v, theta)
  as.numeric(solve(diag(length(v)) + step * S, v))
}

# The parent of a separable penalty, with any fixed() wrapper unwrapped, so
# that the closed forms below can be recognized by the family underneath.
#' The Family Name Under a Separable Penalty's Wrappers
#'
#' @description
#' Returns the bare name of the distribution on which a separable penalty is
#' built. The function removes a `fixed()` wrapper and any bracketed parameter
#' list, so that the closed-form branches can be recognized by the family
#' underneath.
#'
#' @details
#' The parent of a separable penalty is usually wrapped: the lasso's parent is
#' `fixed laplace2 [mu=0]`, and a separable Gaussian prior built as
#' `distrib_penalty(fixed(gaussian1_distrib(), mu = 0), ...)` has the parent
#' `fixed gaussian1 [mu=0]`. A match on the raw `distrib_name` would miss
#' both. Removing the leading `fixed ` and everything from the first bracket
#' leaves `laplace2` and `gaussian1`.
#'
#' The name is used only to reach a closed form. Whether the parent is centered
#' at zero is tested separately, by evaluating the penalty's gradient at the
#' origin, because the family name does not carry the location.
#'
#' @param pen A [DistribPenalty()] object.
#'
#' @return A single string: `"gaussian1"`, `"laplace2"`, `"laplace"` or
#'   `"enet"` for the families with a closed form (the Laplace in two
#'   parametrizations), and the name of the parent's family otherwise.
#'
#' @seealso [penalty_prox.DistribPenalty()], [penalty_prox_spec()]
#'
#' @keywords internal
.prox_family <- function(pen) {
  d <- pen@parent
  nm <- d@distrib_name
  sub("^fixed ", "", sub(" \\[.*$", "", nm))
}

#' The Value of a Parent's Parameter, Free or Held
#'
#' @description
#' Returns the parent distribution's parameter `name`: from the aligned
#' hyperparameter list when the penalty carries it as a free hyperparameter,
#' and from the parent's `fixed_params` when [distributions7::fixed()] holds it
#' at a value.
#'
#' @details
#' The closed forms need the value of the parameter, wherever it is stored. A
#' penalty whose scale is held by [distributions7::fixed()] is the same penalty
#' as one whose scale is free, evaluated at that scale, and
#' [penalty_value()], [penalty_gradient()] and [penalty_hessian()] return the
#' same results for the two. [penalty_prox()] therefore reads the value from
#' `theta` when it is a free hyperparameter and from the parent's
#' `fixed_params` otherwise.
#'
#' @param pen A [DistribPenalty()] object.
#' @param theta The aligned hyperparameter list, as [align_ptheta()] returns it.
#' @param name The parameter's name, a single string.
#'
#' @return A single number.
#'
#' @section Errors:
#' The message names the penalty and the parameter when the parameter is
#' neither a free hyperparameter nor held by [distributions7::fixed()].
#'
#' @seealso [penalty_prox.DistribPenalty()], [penalty_prox_spec()]
#'
#' @keywords internal
.prox_param <- function(pen, theta, name) {
  if (name %in% pen@params) return(theta[[name]])
  fx <- pen@parent@fixed_params
  if (!is.null(fx) && name %in% names(fx)) return(fx[[name]])
  stop(sprintf(paste0(
    "the parent of penalty '%s' carries no parameter '%s', neither as a free\n",
    "  hyperparameter nor held at a value by distributions7::fixed()."),
    pen@penalty_name, name), call. = FALSE)
}

#' @title Proximal Operator of a Separable Penalty
#' @name penalty_prox.DistribPenalty
#'
#' @description
#' Computes closed forms for the Gaussian, the Laplace (in both
#' parametrizations) and the elastic net, and for every other parent with a
#' differentiable log-density it computes the stationary point of smallest
#' objective, coordinate by coordinate. A diagonal map is handled by rescaling,
#' and any other map is rejected.
#'
#' @details
#' # The closed forms
#'
#' | parent | operator |
#' |---|---|
#' | `gaussian1` | \eqn{(v - t g_0)/(1 + t/\sigma^2)}, with \eqn{g_0} the gradient at the origin |
#' | `laplace`, `laplace2` | \eqn{\mathrm{sign}(v)(\lvert v\rvert - t\lambda)_{+}} |
#' | `enet` | the soft threshold at \eqn{t\lambda\alpha}, then division by \eqn{1 + t\lambda(1-\alpha)} |
#'
#' For the `laplace` parametrization \eqn{\lambda = 1/\sigma}. The two steps of
#' the elastic net follow from one stationary condition:
#' \eqn{(\beta - v)/t + a\,\mathrm{sign}(\beta) + c\beta = 0}, with
#' \eqn{a = \lambda\alpha} and \eqn{c = \lambda(1-\alpha)}, separates into the
#' threshold of the Laplace part followed by the shrinkage of the Gaussian
#' part.
#'
#' # The location of the parent
#'
#' The Laplace and elastic-net forms are written for a parent centered at zero,
#' and the family name does not carry the location. The gradient of the
#' penalty at the origin is therefore evaluated, and the call signals an error
#' when its largest absolute entry exceeds \eqn{10^{-8}\max(1, \lambda)}. The
#' Gaussian form has no such restriction, because the location enters the
#' stationary condition linearly and is carried through as \eqn{g_0}.
#'
#' # Parents without a closed form
#'
#' Any other parent is solved coordinate by coordinate from
#' \eqn{(\beta - v)/t = \ell^{(y)}(\beta)}, whose right-hand side is
#' [distributions7::distrib_grad_y()]. The condition can have more than one
#' root when the density is not log-concave, so every change of sign on the
#' segment between \eqn{v} and the mode of the parent is refined and the root
#' of smallest objective is returned, as [.prox_root()] describes. That root
#' is the global minimizer when the parent is unimodal, and a stationary point
#' that need not be the global minimizer when the parent has more than one
#' mode. A parent that declares a kink and is not one of the closed-form
#' families is rejected, because the stationary condition does not hold at a
#' kink and the minimizer of such a penalty can lie there.
#'
#' @param pen A [DistribPenalty()] object.
#' @param v A numeric vector of length `pen@n_coef`.
#' @param step The step length \eqn{t}, a single positive finite number.
#' @param theta A named list of hyperparameter values.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the same length as `v`.
#'
#' @section Errors:
#' A parent read in blocks, a map that is not a \pkg{Matrix} diagonal object
#' with nonzero entries, a Laplace or elastic-net parent that is not centered
#' at zero, a parent that declares a kink and has no closed form, and a
#' stationary condition with no change of sign on the segment or in the
#' widened bracket each signal an error with a message giving the reason; all
#' but the off-center parent name the penalty.
#'
#' @examples
#' v <- c(2, 0.3, -1.4)
#'
#' # The soft threshold sets everything within t * lambda of zero to zero.
#' penalty_prox(lasso_penalty(n_coef = 3), v, 1, list(lambda = 1))
#'
#' # The elastic net thresholds first and then shrinks what survives.
#' penalty_prox(elasticnet_penalty(n_coef = 3), v, 1,
#'              list(lambda = 1, alpha = 0.6))
#' sign(v) * pmax(abs(v) - 1 * 1 * 0.6, 0) / (1 + 1 * 1 * (1 - 0.6))
#'
#' # A Student t prior has no closed form, and the operator is found from the
#' # stationary condition. With sigma = 1, nu = 4 and step 1 it maps 2 to 1,
#' # as the ridge with lambda = 1 does; the two differ for coordinates that
#' # are large relative to sigma, such as 10.
#' v2 <- c(2, 0.3, 10)
#' penalty_prox(heavy_penalty(n_coef = 3), v2, 1, list(sigma = 1, nu = 4))
#' penalty_prox(ridge_penalty(n_coef = 3), v2, 1, list(lambda = 1))
#'
#' @seealso [penalty_prox()] for the generic, [distrib_penalty()] for the
#'   branch, [distributions7::distrib_grad_y()] for the root's right-hand side.
#' @keywords internal
S7::method(penalty_prox, DistribPenalty) <- function(pen, v, step, theta, ...) {
  if (pen@block > 1L) {
    stop(sprintf(paste0(
      "'%s' is read in blocks of %d, and the proximal operator acts one\n",
      "  coordinate at a time: the coordinates of a block do not separate.\n",
      "  has_prox() returns FALSE for it."),
      pen@penalty_name, pen@block), call. = FALSE)
  }
  d <- .prox_scaling(pen)
  if (!is.null(d)) {
    return(S7::method(penalty_prox, DistribPenalty)(
      .undiag(pen), d * v, step * d^2, theta, ...) / d)
  }
  fam <- .prox_family(pen)

  # The closed forms below hold for a parent centered where the quadratic
  # pull is, and a family name does not say where that is. The gradient at
  # the origin does: it is the location term of the stationary condition, so
  # the Gaussian carries it through, and the Laplace, whose location cannot
  # be recovered from one evaluation, is required to sit at zero.
  g0 <- penalty_gradient(pen, rep(0, length(v)), theta)

  if (identical(fam, "gaussian1")) {
    s <- .prox_param(pen, theta, "sigma")
    return((v - step * g0) / (1 + step / s^2))
  }
  if (identical(fam, "laplace2") || identical(fam, "laplace")) {
    lam <- if (identical(fam, "laplace2")) {
      .prox_param(pen, theta, "lambda")
    } else {
      1 / .prox_param(pen, theta, "sigma")
    }
    if (max(abs(g0)) > 1e-8 * max(1, lam)) {
      stop(paste("the Laplace proximal operator is written for a parent",
                 "centered at zero; this one is not."), call. = FALSE)
    }
    return(sign(v) * pmax(abs(v) - step * lam, 0))
  }
  if (identical(fam, "enet")) {
    lam <- .prox_param(pen, theta, "lambda")
    al <- .prox_param(pen, theta, "alpha")
    if (max(abs(g0)) > 1e-8 * max(1, lam)) {
      stop(paste("the elastic-net proximal operator is written for a parent",
                 "centered at zero; this one is not."), call. = FALSE)
    }
    # the stationary condition (b - v)/t + a sgn(b) + c b = 0 separates
    # into the soft threshold of the Laplace part followed by the
    # shrinkage of the Gaussian one
    sft <- sign(v) * pmax(abs(v) - step * lam * al, 0)
    return(sft / (1 + step * lam * (1 - al)))
  }

  if (length(pen@kinks)) {
    stop(sprintf(paste0(
      "'%s' declares a kink and is not one of the closed-form families, so\n",
      "  its proximal operator cannot be found from the stationary condition."),
      pen@penalty_name), call. = FALSE)
  }

  .prox_root(pen, v, step, theta)
}

#' The Proximal Operator of a Separable Penalty, by Its Stationary Points
#'
#' @description
#' Returns, for every coordinate, the stationary point of smallest objective
#' \eqn{(\beta - v_j)^2/(2t) - \log f(\beta)}, where \eqn{f} is the parent's
#' density, among the stationary points found on the segment between
#' \eqn{v_j} and the mode of the parent. When the parent is unimodal, this
#' point is the global minimizer. The function serves the separable penalties
#' that have no closed-form operator, such as the heavy-tailed Student t
#' prior.
#'
#' @details
#' The stationary points are the roots of
#' \eqn{h(\beta) = (\beta - v)/t - \ell^{(y)}(\beta)}, whose second term is
#' [distributions7::distrib_grad_y()]. For a density that is not log-concave
#' \eqn{h} can have more than one root (up to three for the Student t, two of
#' them local minima of the objective), and the first root that a bracket
#' encloses need not be the global minimizer.
#'
#' For a unimodal parent with mode \eqn{m}, the penalty is non-increasing up
#' to \eqn{m} and non-decreasing after it, so the global minimizer lies on the
#' segment between \eqn{v} and \eqn{m}. Moving a point beyond \eqn{v} to
#' \eqn{v}, or a point beyond \eqn{m} to \eqn{m}, lowers the quadratic term
#' and does not raise the penalty. The mode is the root of
#' \eqn{\ell^{(y)}}, found once per call because the parent and \eqn{\theta}
#' are the same for every coordinate: when \eqn{\ell^{(y)}(0)} is exactly
#' zero the mode is taken as zero, and otherwise the root is bracketed by
#' doubling a step from zero in the direction in which \eqn{\ell^{(y)}}
#' points, up to \eqn{2^{60}}, and refined by `uniroot()` to full relative
#' precision. When the search finds no change of sign, or \eqn{\ell^{(y)}} is
#' not finite, the mode is taken as zero and the segment is the one between
#' zero and \eqn{v}. The search moves away from zero only while
#' \eqn{\ell^{(y)}} keeps its sign, so the root that it finds is a local mode,
#' which is the mode when the parent is unimodal.
#'
#' The function scans the segment on `n_grid` equally spaced points that
#' include both ends, and evaluates \eqn{h} for every coordinate in one call.
#' It refines by `uniroot()` every interval over which \eqn{h} changes sign,
#' keeps as it is a grid point where \eqn{h} is exactly zero, and returns the
#' candidate of smallest objective, computing the log-density with
#' [distributions7::distrib_pdf()].
#'
#' When the segment shows no change of sign, the bracket starts one unit
#' beyond \eqn{\min(v, m)} and \eqn{\max(v, m)}, each end that does not yet
#' have the required sign is moved out by the current width of the bracket up
#' to sixty times, and the single root found in the final bracket is
#' returned. For a parent with more than one mode, the result in either case
#' is a stationary point that need not be the global minimizer, because the
#' scan covers only the segment between \eqn{v} and one local mode.
#'
#' The tolerance passed to `uniroot()` is `eps^0.75` times the end of the
#' bracket nearer to zero, or the smallest positive double when the bracket
#' has zero as an end or contains it. `uniroot()` adds `2 eps` times the
#' current iterate to that tolerance, so the root is found to a relative
#' accuracy whatever its magnitude.
#'
#' @param pen A [DistribPenalty()] object with no map and a block width of one.
#' @param v A numeric vector, the points towards which the operator pulls.
#' @param step The step length \eqn{t}, a single positive finite number.
#' @param theta The aligned hyperparameter list.
#' @param n_grid The number of grid points on each segment, including both
#'   ends.
#'
#' @return A numeric vector of the same length as `v`.
#'
#' @return `.prox_mode()` returns a single number, the mode of the parent or
#'   zero.
#'
#' @section Errors:
#' The message names the penalty where the stationary condition changes sign
#' neither on the segment nor in the widened bracket.
#'
#' @seealso [penalty_prox.DistribPenalty()]
#'
#' @keywords internal
.prox_root <- function(pen, v, step, theta, n_grid = 64L) {
  ly <- function(b) {
    as.numeric(distributions7::distrib_grad_y(pen@parent, b, theta)) + 0 * b
  }
  q <- length(v)
  # the global minimizer lies between v_j and the parent's mode m, the same
  # for every coordinate; m = 0 for a parent centered at zero
  m <- .prox_mode(ly)
  # column j of G runs over the segment from m to v_j
  G <- outer(seq(0, 1, length.out = n_grid), v - m) + m
  H <- (G - rep(v, each = n_grid)) / step -
    matrix(ly(as.numeric(G)), nrow = n_grid, ncol = q)
  vapply(seq_len(q), function(j) {
    vj <- v[[j]]
    h <- function(b) (b - vj) / step - ly(b)
    g <- G[, j]
    hj <- H[, j]
    cand <- unique(g[which(hj == 0)])
    for (k in which(hj[-n_grid] * hj[-1L] < 0)) {
      cand <- c(cand, .prox_uniroot(h, g[[k]], g[[k + 1L]],
                                    hj[[k]], hj[[k + 1L]]))
    }
    if (!length(cand)) return(.prox_root_widened(pen, h, vj, m))
    if (length(cand) == 1L) return(cand)
    obj <- (cand - vj)^2 / (2 * step) -
      as.numeric(distributions7::distrib_pdf(pen@parent, cand, theta,
                                             log = TRUE))
    cand[[which.min(obj)]]
  }, numeric(1))
}

#' @rdname dot-prox_root
#' @param ly The function \eqn{\ell^{(y)}} of the parent, vectorized.
#' @keywords internal
.prox_mode <- function(ly) {
  # the root of l_y bracketed by doubling from 0 in the direction in which
  # l_y points; when l_y does not change sign within 2^60, or is not finite,
  # the mode is taken as 0 and the scan runs from 0 to v as for a centered
  # parent
  f0 <- ly(0)
  if (!is.finite(f0) || f0 == 0) return(0)
  a <- 0
  fa <- f0
  for (k in 0:60) {
    b <- sign(f0) * 2^k
    fb <- ly(b)
    if (!is.finite(fb)) return(0)
    if (fb == 0) return(b)
    if (sign(fb) != sign(f0)) {
      # to full precision: the root of h next to the mode can lie much
      # closer to it than one scale, and it has to stay on the segment
      lo <- min(a, b)
      hi <- max(a, b)
      return(stats::uniroot(ly, c(lo, hi),
                            f.lower = if (lo == a) fa else fb,
                            f.upper = if (hi == a) fa else fb,
                            tol = .Machine$double.xmin, maxiter = 2000L)$root)
    }
    a <- b
    fa <- fb
  }
  0
}

#' @rdname dot-prox_root
#' @param h The function \eqn{h} of one coordinate.
#' @param vi That coordinate's point.
#' @param m The parent's mode, as `.prox_mode()` returns it.
#' @keywords internal
.prox_root_widened <- function(pen, h, vi, m = 0) {
  lo <- min(vi, m) - 1
  hi <- max(vi, m) + 1
  for (k in 1:60) {
    if (h(lo) <= 0 && h(hi) >= 0) break
    if (h(lo) > 0) lo <- lo - (hi - lo)
    if (h(hi) < 0) hi <- hi + (hi - lo)
  }
  hl <- h(lo)
  hh <- h(hi)
  if (hl > 0 || hh < 0) {
    stop(sprintf(paste0(
      "the proximal operator of '%s' could not be found: the stationary\n",
      "  condition changes sign neither between %s and v nor in a bracket\n",
      "  widened sixty times."), pen@penalty_name,
      if (m == 0) "0" else "the mode of the parent"), call. = FALSE)
  }
  if (hl == 0) return(lo)
  if (hh == 0) return(hi)
  .prox_uniroot(h, lo, hi, hl, hh)
}

#' @rdname dot-prox_root
#' @param lo,hi The ends of a bracket, in either order.
#' @param flo,fhi The values of `h` there, of opposite signs.
#' @keywords internal
.prox_uniroot <- function(h, lo, hi, flo, fhi) {
  if (lo > hi) {
    b <- c(hi, lo)
    f <- c(fhi, flo)
    lo <- b[[1L]]
    hi <- b[[2L]]
    flo <- f[[1L]]
    fhi <- f[[2L]]
  }
  # uniroot() stops once the bracket is within 2 eps |b| + tol/2, so a tol
  # proportional to the end nearer zero is relative to the root
  m <- if (lo > 0 || hi < 0) min(abs(lo), abs(hi)) else 0
  stats::uniroot(h, c(lo, hi), f.lower = flo, f.upper = fhi,
                 tol = max(.Machine$double.eps^0.75 * m,
                           .Machine$double.xmin))$root
}

#' @title Proximal Operator of SCAD
#' @name penalty_prox.ScadPenalty
#'
#' @description
#' Computes the closed piecewise operator of Fan and Li (2001): a soft
#' threshold near zero, a rescaled threshold over the tapering region, and the
#' identity beyond \eqn{a\lambda}, so a large coefficient passes through
#' untouched.
#'
#' @details
#' With \eqn{u = \lvert v\rvert} and \eqn{s} its sign,
#'
#' \deqn{\mathrm{prox}(v) =
#'   \begin{cases}
#'     s\,(u - t\lambda)_{+} & u \le (1+t)\lambda \\
#'     \dfrac{s\,(u - t a \lambda/(a-1))}{1 - t/(a-1)}
#'       & (1+t)\lambda < u \le a\lambda \\
#'     v & u > a\lambda.
#'   \end{cases}}
#'
#' The middle piece is computed in the equivalent form
#' \eqn{s\,(u - t(a\lambda - u)/(a - 1 - t))}, in which \eqn{a\lambda - u} and
#' \eqn{a - 1 - t} are differences of exact operands, so that the value stays
#' accurate as \eqn{t} approaches \eqn{a - 1}.
#'
#' The step enters through the division by \eqn{1 - t/(a-1)}: the curvature of
#' SCAD is \eqn{-1/(a-1)} in the middle piece, so the subproblem is convex only
#' while \eqn{t < a-1}. At or beyond that bound the operator is set-valued, and
#' the call signals an error that prints \eqn{a - 1} and the largest scaled
#' step, defined below.
#'
#' The formula is written for unit curvature. With the curvature \eqn{c_j} of
#' [scad_penalty()] (argument `curv`), the penalty in coordinate \eqn{j} is
#' \eqn{c_j\rho(\cdot;\ \lambda/c_j)}, and the operator is the formula above
#' with \eqn{t c_j} in place of \eqn{t} and \eqn{\lambda/c_j} in place of
#' \eqn{\lambda}. Under a diagonal map the step is \eqn{t d_j^2}, so the
#' condition becomes \eqn{t c_j d_j^2 < a - 1} in every coordinate, which is
#' stricter than \eqn{t < a - 1} where \eqn{c_j d_j^2 > 1} and looser where
#' \eqn{c_j d_j^2 < 1}.
#'
#' @param pen A [ScadPenalty()] object.
#' @param v A numeric vector of length `pen@n_coef`.
#' @param step The step length \eqn{t}, a single positive finite number with
#'   \eqn{t c_j < a-1} in every coordinate (\eqn{t < a - 1} at unit
#'   curvature).
#' @param theta A named list holding `lambda` and `a`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the same length as `v`.
#'
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. *Journal of the American Statistical
#' Association* **96**, 1348-1360.
#'
#' @examples
#' th <- list(lambda = 1, a = 3.7)
#'
#' # The three points lie one in each region: the first is thresholded, the
#' # second rescaled and the third untouched.
#' penalty_prox(scad_penalty(n_coef = 3), c(1.5, 3, 6), 1, th)
#'
#' # Past a * lambda the operator is the identity, so a large coefficient is
#' # not shrunk at all. The lasso shrinks the same point by t * lambda.
#' penalty_prox(scad_penalty(n_coef = 1), 6, 1, th)
#' penalty_prox(lasso_penalty(n_coef = 1), 6, 1, list(lambda = 1))
#'
#' # At the convexity bound the operator is set-valued, and the call signals
#' # an error.
#' try(penalty_prox(scad_penalty(n_coef = 1), 2, 2.7, th))
#'
#' @seealso [penalty_prox()] for the generic, [scad_penalty()] for the penalty,
#'   [penalty_prox.McpPenalty()] for the sibling family,
#'   [penalty_prox_spec()] for the same operator as a table.
#' @keywords internal
S7::method(penalty_prox, ScadPenalty) <- function(pen, v, step, theta, ...) {
  d <- .prox_scaling(pen)
  if (!is.null(d)) {
    return(S7::method(penalty_prox, ScadPenalty)(
      .undiag(pen), d * v, step * d^2, theta, ...) / d)
  }
  cv <- curv_of(pen, length(v))
  lam <- theta$lambda / cv
  a <- theta$a
  step <- step * cv
  if (any(step >= a - 1)) {
    stop(sprintf(paste0(
      "the SCAD proximal operator needs step < a - 1 (%g here), where the\n",
      "  step is multiplied by the curvature and, under a diagonal map, by\n",
      "  the squared diagonal entry; the largest scaled step is %g. Beyond\n",
      "  that bound the subproblem is not convex and the operator is\n",
      "  set-valued."),
      a - 1, max(step)),
      call. = FALSE)
  }
  s <- sign(v)
  u <- abs(v)
  # the middle piece as u less a correction: near t = a - 1 both a lam - u
  # and a - 1 - t are exact differences, where 1 - t/(a - 1) rounds first
  ifelse(u <= (1 + step) * lam,
         s * pmax(u - step * lam, 0),
         ifelse(u <= a * lam,
                s * (u - step * (a * lam - u) / (a - 1 - step)),
                v))
}

#' @title Proximal Operator of MCP
#' @name penalty_prox.McpPenalty
#'
#' @description
#' Computes the closed piecewise operator of Zhang (2010): a rescaled soft
#' threshold below \eqn{\gamma\lambda} and the identity beyond it, so a large
#' coefficient passes through untouched.
#'
#' @details
#' With \eqn{u = \lvert v\rvert} and \eqn{s} its sign,
#'
#' \deqn{\mathrm{prox}(v) =
#'   \begin{cases}
#'     \dfrac{s\,(u - t\lambda)_{+}}{1 - t/\gamma}
#'       & u \le \gamma\lambda \\
#'     v & u > \gamma\lambda.
#'   \end{cases}}
#'
#' MCP has two pieces where SCAD has three, because its taper starts at the
#' origin. The first piece is computed as
#' \eqn{s\,(u - t\lambda)_{+}\,\gamma/(\gamma - t)}, in which
#' \eqn{\gamma - t} is a difference of exact operands, so that the value stays
#' accurate as \eqn{t} approaches \eqn{\gamma}.
#'
#' The step enters through the division by \eqn{1 - t/\gamma}: the curvature
#' is the constant \eqn{-1/\gamma} inside the taper, so the subproblem is
#' convex only while \eqn{t < \gamma}. At or beyond that bound the operator is
#' set-valued, and the call signals an error that prints \eqn{\gamma} and the
#' largest scaled step, defined below.
#'
#' The formula is written for unit curvature. With the curvature \eqn{c_j} of
#' [mcp_penalty()] (argument `curv`), the operator is the formula above with
#' \eqn{t c_j} in place of \eqn{t} and \eqn{\lambda/c_j} in place of
#' \eqn{\lambda}. Under a diagonal map the step is \eqn{t d_j^2}, so the
#' condition becomes \eqn{t c_j d_j^2 < \gamma} in every coordinate, which is
#' stricter than \eqn{t < \gamma} where \eqn{c_j d_j^2 > 1} and looser where
#' \eqn{c_j d_j^2 < 1}.
#'
#' @param pen An [McpPenalty()] object.
#' @param v A numeric vector of length `pen@n_coef`.
#' @param step The step length \eqn{t}, a single positive finite number with
#'   \eqn{t c_j < \gamma} in every coordinate (\eqn{t < \gamma} at unit
#'   curvature).
#' @param theta A named list holding `lambda` and `gamma`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the same length as `v`.
#'
#' @references
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. *Annals of Statistics* **38**, 894-942.
#'
#' @examples
#' th <- list(lambda = 1, gamma = 3)
#'
#' # Inside the taper the soft-thresholded value is divided by 1 - t/gamma;
#' # beyond gamma * lambda the point is returned unchanged.
#' penalty_prox(mcp_penalty(n_coef = 3), c(1.5, 2.5, 6), 1, th)
#'
#' # At the same lambda, MCP with gamma = 3 shrinks this point less than SCAD
#' # with a = 3.7.
#' penalty_prox(mcp_penalty(n_coef = 1), 2.5, 1, th)
#' penalty_prox(scad_penalty(n_coef = 1), 2.5, 1, list(lambda = 1, a = 3.7))
#'
#' # At the convexity bound the operator is set-valued, and the call signals
#' # an error.
#' try(penalty_prox(mcp_penalty(n_coef = 1), 2, 3, th))
#'
#' @seealso [penalty_prox()] for the generic, [mcp_penalty()] for the penalty,
#'   [penalty_prox.ScadPenalty()] for the sibling family,
#'   [penalty_prox_spec()] for the same operator as a table.
#' @keywords internal
S7::method(penalty_prox, McpPenalty) <- function(pen, v, step, theta, ...) {
  d <- .prox_scaling(pen)
  if (!is.null(d)) {
    return(S7::method(penalty_prox, McpPenalty)(
      .undiag(pen), d * v, step * d^2, theta, ...) / d)
  }
  cv <- curv_of(pen, length(v))
  lam <- theta$lambda / cv
  gam <- theta$gamma
  step <- step * cv
  if (any(step >= gam)) {
    stop(sprintf(paste0(
      "the MCP proximal operator needs step < gamma (%g here), where the\n",
      "  step is multiplied by the curvature and, under a diagonal map, by\n",
      "  the squared diagonal entry; the largest scaled step is %g. Beyond\n",
      "  that bound the subproblem is not convex and the operator is\n",
      "  set-valued."),
      gam, max(step)),
      call. = FALSE)
  }
  s <- sign(v)
  u <- abs(v)
  # gam - t is an exact difference near t = gam, where 1 - t/gam rounds first
  ifelse(u <= gam * lam,
         s * pmax(u - step * lam, 0) * gam / (gam - step),
         v)
}
