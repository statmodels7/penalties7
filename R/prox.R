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
#' A proximal gradient method applies it after each gradient step. It is the
#' operation that lets a penalty with a kink be minimized without
#' differencing it, and the one that sets a coefficient exactly to zero where a
#' gradient method only makes it small.
#'
#' @details
#' # The four routes
#'
#' | branch | operator |
#' |---|---|
#' | quadratic, structured | one linear solve, \eqn{(I + tS)^{-1}v} |
#' | Gaussian prior | the shrinkage \eqn{v/(1 + t/\sigma^2)} |
#' | Laplace prior | the soft threshold \eqn{\mathrm{sign}(v)(\lvert v\rvert - t\lambda)_{+}} |
#' | elastic net | the soft threshold followed by the shrinkage |
#' | SCAD, MCP | the closed piecewise operators of their papers |
#' | any other separable parent | a coordinatewise root, found by bisection |
#'
#' Checked against a direct minimization of the defining objective,
#' coordinate by coordinate on a three-vector at a step of 1, all eight shipped
#' penalties agree to `2.7e-08` or better, which is the accuracy of the
#' minimizer used as the reference and not of the operator.
#'
#' # The root, where there is no formula
#'
#' Any separable penalty whose parent has a differentiable log-density is
#' solved from the stationary condition
#'
#' \deqn{\frac{\beta - v}{t} = \ell^{(y)}(\beta),}
#'
#' with the right-hand side [distributions7::distrib_grad_y()], closed form for
#' every continuous family. Log-concavity of the density makes the left side
#' minus the right side strictly increasing, so the root is unique; the
#' bracket is widened up to sixty times and then `uniroot()` finds it to
#' `eps^0.75`. A parent that declares a kink and is not one of the closed-form
#' families is rejected instead, a smooth root being the wrong instrument
#' there.
#'
#' # Under a map
#'
#' The quadratic and structured branches take **any** map, the objective
#' staying quadratic. The separable branches need the identity map or a
#' diagonal one, under which the operator is the unmapped one read at
#' \eqn{dv} with the step \eqn{td^2}, divided back by \eqn{d}: that is what
#' standardization comes to. A general \eqn{D} mixes coordinates and makes this
#' the generalized-lasso problem, which needs an algorithm of its own, so it is
#' rejected by name.
#'
#' # What rejects
#'
#' - a penalty with no operator at all, [additive_penalty()] being the one that
#'   ships, with a message naming [has_prox()];
#' - a separable penalty under a map that is not diagonal;
#' - a separable penalty whose parent is read in blocks, whose coordinates do
#'   not separate;
#' - a Laplace or elastic-net parent not centered at zero, the closed form
#'   being written for one that is;
#' - SCAD at `step >= a - 1` and MCP at `step >= gamma`, where the subproblem
#'   is not convex and the operator is set-valued;
#' - a `step` that is not a single positive number, or a `v` whose length is
#'   not `pen@n_coef`.
#'
#' @param pen A [penalty()] object. Ask [has_prox()] first if you do not know
#'   the branch.
#' @param v A numeric vector of length `pen@n_coef`, the point pulled towards.
#'   Coerced with `as.numeric()`.
#' @param step The step length \eqn{t}, a single positive number. A vector is
#'   rejected; [penalty_prox_spec()] is the interface that takes one per
#'   coefficient.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, holding every name in `pen@params`.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A numeric vector of the same length as `v`.
#'
#' @examples
#' v <- c(2, 0.3, -1.4)
#'
#' # The lasso: a soft threshold, so the middle coordinate becomes exactly 0.
#' penalty_prox(lasso_penalty(n_coef = 3), v, 1, list(lambda = 1))
#'
#' # A ridge: a shrinkage towards zero, and nothing reaches it.
#' penalty_prox(ridge_penalty(n_coef = 3), v, 1, list(lambda = 1.5))
#' v / (1 + 1 * 1.5)
#'
#' # The elastic net is the one followed by the other.
#' penalty_prox(elasticnet_penalty(n_coef = 3), v, 1,
#'              list(lambda = 1, alpha = 0.6))
#'
#' # A heavy-tailed prior has no closed form, so this is a bracketed root.
#' # It leaves the large coordinate nearly alone, which is the point of it.
#' penalty_prox(heavy_penalty(n_coef = 3), v, 1, list(sigma = 1, nu = 4))
#'
#' # Every one of them minimizes the defining objective. Checking the lasso's
#' # first coordinate directly:
#' rho <- lasso_penalty(n_coef = 3)
#' obj <- function(b1) {
#'   x <- v; x[1] <- b1
#'   0.5 * sum((x - v)^2) / 1 + penalty_value(rho, x, list(lambda = 1))
#' }
#' stats::optimize(obj, c(-5, 5), tol = 1e-12)$minimum
#'
#' # A branch with no operator says so rather than approximating.
#' try(penalty_prox(additive_penalty(list(diag(3), diag(c(1, 1, 0)))),
#'                  v, 1, list(lambda1 = 1, lambda2 = 1)))
#'
#' @seealso [has_prox()] to ask before calling, [penalty_prox_spec()] for the
#'   same operator as a table a compiled loop can read, [penalty_kinks()] for
#'   where a gradient method would fail, [penalty_gradient()] for what to use
#'   where the penalty is smooth.
#'
#' @export
penalty_prox <- S7::new_generic("penalty_prox", "pen",
  function(pen, v, step, theta, ...) {
    theta <- align_ptheta(pen, theta)
    v <- as.numeric(v)
    if (length(step) != 1L || is.na(step) || step <= 0) {
      stop("'step' must be a single positive number.", call. = FALSE)
    }
    if (length(v) != pen@n_coef) {
      stop(sprintf("'v' must have length %d.", pen@n_coef), call. = FALSE)
    }
    S7::S7_dispatch()
  })

S7::method(penalty_prox, penalty) <- function(pen, v, step, theta, ...) {
  stop(sprintf(paste0(
    "'%s' has no proximal operator.\n",
    "  A penalty supplies one when it is quadratic, or separable under the\n",
    "  identity map. Use has_prox() to ask before calling."),
    pen@penalty_name), call. = FALSE)
}

#' @title Does a Penalty Supply a Proximal Operator?
#'
#' @description
#' `TRUE` when [penalty_prox()] can be evaluated for this penalty, so that a
#' caller may route a block to a proximal method without provoking an error.
#' `FALSE` when the penalty has none, and then the block belongs to a smooth
#' method or to a scheme of its own.
#'
#' @details
#' # What it tests
#'
#' The predicate asks four things in turn: whether the class registers a
#' [penalty_prox()] method of its own rather than inheriting the base class's
#' refusal; whether the penalty is quadratic or structured, in which case the
#' operator is a linear solve at any map; whether a separable parent is read in
#' blocks, whose coordinates do not separate; and otherwise whether the map is
#' absent or diagonal.
#'
#' Measured over the shipped branches:
#'
#' | penalty | `has_prox()` |
#' |---|---|
#' | quadratic, at any map | `TRUE` |
#' | structured | `TRUE` |
#' | ridge, lasso, elastic net, heavy-tailed | `TRUE` |
#' | lasso under a **diagonal** map | `TRUE` |
#' | SCAD, MCP | `TRUE` |
#' | **additive** | `FALSE` |
#' | lasso under a **general** map | `FALSE` |
#' | a separable penalty read in blocks | `FALSE` |
#'
#' # What `TRUE` does not promise
#'
#' It says the penalty has an operator, not that every step reaches it. SCAD
#' and MCP answer `TRUE` and still reject a step at or beyond \eqn{a - 1} and
#' \eqn{\gamma}, where the subproblem is not convex and the operator is
#' set-valued. A Laplace or elastic-net parent that is not centered at zero
#' answers `TRUE` and rejects at the call. The predicate is about the penalty;
#' those two conditions are about the step and the parent.
#'
#' @param pen A [penalty()] object of any branch. Anything else is rejected.
#'
#' @return A single logical.
#'
#' @examples
#' c(lasso = has_prox(lasso_penalty()),
#'   scad  = has_prox(scad_penalty()),
#'   ridge = has_prox(ridge_penalty()))
#'
#' # The additive branch is the one that ships without an operator.
#' has_prox(additive_penalty(list(diag(3), diag(c(1, 1, 0)))))
#'
#' # A map decides it for a separable penalty: diagonal yes, general no.
#' c(diagonal = has_prox(lasso_penalty(map = Matrix::Diagonal(x = c(1, 2, 3)))),
#'   general  = has_prox(lasso_penalty(map = rbind(c(1, -1, 0), c(0, 1, -1)))))
#'
#' # TRUE does not mean every step works: SCAD needs step < a - 1.
#' has_prox(scad_penalty(n_coef = 1))
#' try(penalty_prox(scad_penalty(n_coef = 1), 2, 2.7, list(lambda = 1, a = 3.7)))
#'
#' @seealso [penalty_prox()] for the operator, [penalty_prox_spec()] for it as
#'   a table, [is_quadratic()] and [is_proper()] for the other two predicates a
#'   consumer routes on.
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
#' Returns the diagonal entries of \eqn{D} where the map is a \pkg{Matrix}
#' diagonal object, and `NULL` where there is no map, where the map is of
#' another class, or where a diagonal entry is zero or missing.
#'
#' @details
#' A diagonal map rescales each coordinate on its own, and a separable penalty
#' under one is still separable. That is what standardization comes to:
#' penalizing a column divided by its own spread is penalizing
#' \eqn{\rho(s_j\beta_j)}, so the scaling never has to touch the design and the
#' sparsity of a block survives it. A general \eqn{D} mixes coordinates and
#' turns the problem into the generalized-lasso one, which is a different
#' algorithm and not a different formula.
#'
#' The map is recognized by its **class** and not by inspecting its entries: a
#' \pkg{Matrix} diagonal object says what it is and costs \eqn{q} numbers,
#' where testing a dense matrix for diagonality would cost \eqn{q^2} and defeat
#' the point. A base matrix that happens to be diagonal therefore answers
#' `NULL`, and a penalty carrying one has no proximal operator.
#'
#' A zero on the diagonal is rejected as well, the transport dividing by it.
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
#' Returns the same penalty with `map` set to `NULL`, which is the object the
#' diagonal-map transport hands to the identity-map formulas.
#'
#' @details
#' The transport applies the diagonal to the point and to the step before
#' calling, and divides the answer back afterwards, so what the formula needs
#' is the penalty without its map. Nothing else on the object changes, and in
#' particular `n_coef` is unaffected, a diagonal map being square.
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

#' The Diagonal a Proximal Operator Rescales By
#'
#' @description
#' Returns the map's diagonal where there is one, `NULL` where there is no map,
#' and an error naming the penalty where the map is not diagonal. The three
#' answers are the three things a proximal method must do: rescale, proceed,
#' or decline.
#'
#' @details
#' The identity behind the rescaling is
#'
#' \deqn{\mathrm{prox}_{t\rho(d\,\cdot)}(v) =
#'   \mathrm{prox}_{t d^2 \rho}(d v)/d,}
#'
#' so a caller with a diagonal map applies \eqn{d} to the point, squares it
#' into the step, calls the identity-map formula and divides back.
#'
#' The error is raised here, at the call, rather than left to a formula that
#' would silently ignore the map. [has_prox()] answers `FALSE` for the same
#' penalty, so a caller that asks first never reaches it.
#'
#' @param pen A [penalty()] object.
#'
#' @return A numeric vector of diagonal entries, or `NULL` when the map is
#'   absent.
#'
#' @section Errors:
#' A message naming the penalty and saying that a non-diagonal map makes this
#' the generalized-lasso problem, where the operator does not split by
#' coordinate.
#'
#' @seealso [map_diagonal()], [.undiag()], [penalty_prox()]
#'
#' @keywords internal
.prox_scaling <- function(pen) {
  if (is.null(pen@map)) return(NULL)
  d <- map_diagonal(pen)
  if (is.null(d)) {
    stop(sprintf(paste0(
      "'%s' carries a linear map that is not diagonal, and its proximal\n",
      "  operator does not split by coordinate. Only the identity map and a\n",
      "  diagonal one have a closed form here."),
      pen@penalty_name), call. = FALSE)
  }
  d
}

#' @title Proximal Operator of a Quadratic or Structured Penalty
#' @name penalty_prox.QuadraticPenalty
#'
#' @description
#' One linear solve, \eqn{(I + tS)^{-1}v}, with \eqn{S} the penalty's Hessian:
#' \eqn{\lambda D'PD} for a quadratic penalty and \eqn{\Omega(\theta)} for a
#' structured one. The two branches share this body, the objective being
#' quadratic in both.
#'
#' @details
#' The subproblem is
#' \eqn{\arg\min_\beta \{\tfrac{1}{2t}\lVert\beta - v\rVert^2 +
#' \tfrac{1}{2}\beta'S\beta\}}, whose stationary condition is
#' \eqn{(\beta - v)/t + S\beta = 0}. The normalizing constant does not depend
#' on \eqn{\beta} and drops out.
#'
#' Because the objective stays quadratic under any linear map, this route needs
#' no restriction on \eqn{D}, where the separable branches need the identity or
#' a diagonal. The solve is dense and costs \eqn{O(q^3)} whatever the storage
#' of \eqn{S}, which for a penalty built with `blocks > 1` is a `dgCMatrix`
#' that is densified by `diag(length(v)) + step * S`.
#'
#' @param pen A [QuadraticPenalty()] or [StructuredPenalty()] object.
#' @param v A numeric vector of length `pen@n_coef`.
#' @param step The step length \eqn{t}, a single positive number.
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
#' # A second-difference penalty pulls towards the straight lines it does not
#' # charge for, rather than towards zero.
#' penalty_prox(quadratic_penalty(crossprod(diff(diag(3)))), v, 1,
#'              list(lambda = 5))
#'
#' # A structured penalty at a zero log-Cholesky free vector is the ridge.
#' s <- structured_penalty(parameters7::log_cholesky(3, role = "precision"))
#' penalty_prox(s, v, 1,
#'              as.list(stats::setNames(rep(0, 6), s@params)))
#'
#' @seealso [penalty_prox()] for the generic and the other branches,
#'   [penalty_hessian()] for the matrix solved against.
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
#' Returns the bare name of the distribution a separable penalty is built on,
#' with a `fixed()` wrapper and any bracketed parameter list stripped off, so
#' that the closed-form branches can be recognized by the family underneath.
#'
#' @details
#' A separable penalty's parent is typically wrapped: the lasso is
#' `fixed laplace2 [mu=0]` and the ridge's separable twin is
#' `fixed gaussian1 [mu=0]`. Matching on the raw `distrib_name` would miss both.
#' Stripping the leading `fixed ` and everything from the first bracket leaves
#' `laplace2` and `gaussian1`.
#'
#' The name is used only to reach a closed form. Whether the parent is centered
#' where the quadratic pull is is asked separately, by evaluating the gradient
#' at the origin, since a family name does not say where a location sits.
#'
#' @param pen A [DistribPenalty()] object.
#'
#' @return A single string: `"gaussian1"`, `"laplace2"`, `"laplace"`, `"enet"`
#'   for the four closed-form families, and whatever the parent is called
#'   otherwise.
#'
#' @seealso [penalty_prox.DistribPenalty()], [penalty_prox_spec()]
#'
#' @keywords internal
.prox_family <- function(pen) {
  d <- pen@parent
  nm <- d@distrib_name
  sub("^fixed ", "", sub(" \\[.*$", "", nm))
}

#' @title Proximal Operator of a Separable Penalty
#' @name penalty_prox.DistribPenalty
#'
#' @description
#' Closed form for the three families that have one, and a coordinatewise root
#' for every other parent with a differentiable log-density. A diagonal map is
#' handled by rescaling; anything else is rejected.
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
#' The elastic net's two steps come out of one stationary condition:
#' \eqn{(\beta - v)/t + a\,\mathrm{sign}(\beta) + c\beta = 0} with
#' \eqn{a = \lambda\alpha} and \eqn{c = \lambda(1-\alpha)} separates into the
#' Laplace part's threshold followed by the Gaussian part's shrinkage.
#'
#' # Where the parent must sit
#'
#' The Laplace and elastic-net forms are written for a parent centered at zero,
#' and a family name does not say where a location is. The gradient at the
#' origin is evaluated and required to vanish, and an off-center parent is
#' rejected. The Gaussian needs no such restriction: its location enters the
#' stationary condition linearly and is carried through as \eqn{g_0}.
#'
#' # The root
#'
#' Any other parent is solved coordinatewise from
#' \eqn{(\beta - v)/t = \ell^{(y)}(\beta)}, the right-hand side being
#' [distributions7::distrib_grad_y()]. For a log-concave density the difference
#' of the two sides is strictly increasing, so the root is unique; the bracket
#' starts one unit either side of \eqn{\min(v, 0)} and \eqn{\max(v, 0)} and is
#' doubled up to sixty times, then `uniroot()` refines to `eps^0.75`. A parent
#' declaring a kink and not among the closed-form families is rejected instead:
#' the root is a smooth instrument and the kink is where the answer is.
#'
#' @param pen A [DistribPenalty()] object.
#' @param v A numeric vector of length `pen@n_coef`.
#' @param step The step length \eqn{t}, a single positive number.
#' @param theta A named list of hyperparameter values.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector of the same length as `v`.
#'
#' @section Errors:
#' A parent read in blocks, a map that is not diagonal, an off-center Laplace
#' or elastic-net parent, a parent with an unrecognized kink, and a bracket
#' that could not be found, each with a message naming the penalty and the
#' reason.
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
#' # A Student t prior has no closed form. The root leaves the large
#' # coordinate nearly where it was, which a Gaussian prior would not.
#' penalty_prox(heavy_penalty(n_coef = 3), v, 1, list(sigma = 1, nu = 4))
#' penalty_prox(ridge_penalty(n_coef = 3), v, 1, list(lambda = 1))
#'
#' @seealso [penalty_prox()] for the generic, [distrib_penalty()] for the
#'   branch, [distributions7::distrib_grad_y()] for the root's right-hand side.
#' @keywords internal
S7::method(penalty_prox, DistribPenalty) <- function(pen, v, step, theta, ...) {
  if (pen@block > 1L) {
    stop(sprintf(paste0(
      "'%s' is read in blocks of %d, and the proximal operator acts one\n",
      "  coordinate at a time: the coordinates of a block do not separate.\n",
      "  has_prox() answers FALSE for it."),
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
    s <- theta[[which(pen@params == "sigma")]]
    return((v - step * g0) / (1 + step / s^2))
  }
  if (identical(fam, "laplace2") || identical(fam, "laplace")) {
    lam <- if (identical(fam, "laplace2")) {
      theta[[which(pen@params == "lambda")]]
    } else {
      1 / theta[[which(pen@params == "sigma")]]
    }
    if (max(abs(g0)) > 1e-8 * max(1, lam)) {
      stop(paste("the Laplace proximal operator is written for a parent",
                 "centered at zero; this one is not."), call. = FALSE)
    }
    return(sign(v) * pmax(abs(v) - step * lam, 0))
  }
  if (identical(fam, "enet")) {
    lam <- theta[[which(pen@params == "lambda")]]
    al <- theta[[which(pen@params == "alpha")]]
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
      "  its proximal operator cannot be found by a smooth root."),
      pen@penalty_name), call. = FALSE)
  }

  # h(b) = (b - v)/t - l_y(b) is strictly increasing for a log-concave
  # density, so one bracketed root per coordinate.
  ly <- function(b) {
    as.numeric(distributions7::distrib_grad_y(pen@parent, b, theta)) + 0 * b
  }
  vapply(v, function(vi) {
    h <- function(b) (b - vi) / step - ly(b)
    lo <- min(vi, 0) - 1
    hi <- max(vi, 0) + 1
    for (k in 1:60) {
      if (h(lo) <= 0 && h(hi) >= 0) break
      if (h(lo) > 0) lo <- lo - (hi - lo)
      if (h(hi) < 0) hi <- hi + (hi - lo)
    }
    if (h(lo) > 0 || h(hi) < 0) {
      stop("the proximal root could not be bracketed; is the density log-concave?",
           call. = FALSE)
    }
    stats::uniroot(h, c(lo, hi), tol = .Machine$double.eps^0.75)$root
  }, numeric(1))
}

#' @title Proximal Operator of SCAD
#' @name penalty_prox.ScadPenalty
#'
#' @description
#' The closed piecewise operator of Fan and Li (2001): a soft threshold near
#' zero, a rescaled threshold over the tapering region, and the identity beyond
#' \eqn{a\lambda}, so a large coefficient passes through untouched.
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
#' The middle piece divides by \eqn{1 - t/(a-1)}, which is where the step
#' matters: SCAD's curvature is \eqn{-1/(a-1)} there, so the subproblem is
#' convex only while \eqn{t < a-1}. At or beyond that the operator is
#' set-valued and the call is rejected with the bound printed. Under a diagonal
#' map the condition tightens to \eqn{t < (a-1)/d_j^2}.
#'
#' @param pen A [ScadPenalty()] object.
#' @param v A numeric vector of length `pen@n_coef`.
#' @param step The step length \eqn{t}, a single positive number below
#'   \eqn{a-1}.
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
#' # One point in each region: thresholded, rescaled, untouched.
#' penalty_prox(scad_penalty(n_coef = 3), c(1.5, 3, 6), 1, th)
#'
#' # Past a * lambda the operator is the identity, so a large coefficient is
#' # not shrunk at all. The lasso shrinks the same point by t * lambda.
#' penalty_prox(scad_penalty(n_coef = 1), 6, 1, th)
#' penalty_prox(lasso_penalty(n_coef = 1), 6, 1, list(lambda = 1))
#'
#' # Beyond the convexity bound the operator is set-valued and rejects.
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
  lam <- theta$lambda
  a <- theta$a
  if (any(step >= a - 1)) {
    stop(sprintf(paste0(
      "the SCAD proximal operator needs step < a - 1 (%g here): beyond that the\n",
      "  subproblem is not convex and the operator is set-valued."), a - 1),
      call. = FALSE)
  }
  s <- sign(v)
  u <- abs(v)
  ifelse(u <= (1 + step) * lam,
         s * pmax(u - step * lam, 0),
         ifelse(u <= a * lam,
                s * (u - step * a * lam / (a - 1)) / (1 - step / (a - 1)),
                v))
}

#' @title Proximal Operator of MCP
#' @name penalty_prox.McpPenalty
#'
#' @description
#' The closed piecewise operator of Zhang (2010): a rescaled soft threshold
#' below \eqn{\gamma\lambda} and the identity beyond it, so a large coefficient
#' passes through untouched.
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
#' MCP has two pieces where SCAD has three, its taper starting at the origin.
#' The division by \eqn{1 - t/\gamma} is where the step matters: the curvature
#' is a constant \eqn{-1/\gamma} inside the taper, so the subproblem is convex
#' only while \eqn{t < \gamma}. At or beyond that the operator is set-valued
#' and the call is rejected with the bound printed. Under a diagonal map the
#' condition tightens to \eqn{t < \gamma/d_j^2}.
#'
#' @param pen An [McpPenalty()] object.
#' @param v A numeric vector of length `pen@n_coef`.
#' @param step The step length \eqn{t}, a single positive number below
#'   \eqn{\gamma}.
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
#' # Inside the taper the threshold is rescaled upward by 1/(1 - t/gamma);
#' # beyond gamma * lambda the point is returned unchanged.
#' penalty_prox(mcp_penalty(n_coef = 3), c(1.5, 2.5, 6), 1, th)
#'
#' # The rescaling lets MCP reach the unpenalized answer sooner than SCAD at
#' # the same lambda.
#' penalty_prox(mcp_penalty(n_coef = 1), 2.5, 1, th)
#' penalty_prox(scad_penalty(n_coef = 1), 2.5, 1, list(lambda = 1, a = 3.7))
#'
#' # Beyond the convexity bound the operator is set-valued and rejects.
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
  lam <- theta$lambda
  gam <- theta$gamma
  if (any(step >= gam)) {
    stop(sprintf(paste0(
      "the MCP proximal operator needs step < gamma (%g here): beyond that the\n",
      "  subproblem is not convex and the operator is set-valued."), gam),
      call. = FALSE)
  }
  s <- sign(v)
  u <- abs(v)
  ifelse(u <= gam * lam,
         s * pmax(u - step * lam, 0) / (1 - step / gam),
         v)
}
