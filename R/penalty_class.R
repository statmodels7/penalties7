#' @title S7 Base Class for Penalties
#'
#' @description
#' The abstract parent of every penalty in this package. A penalty is a scalar
#' function of the coefficients, written \eqn{\rho(D\beta; \theta)}: a linear
#' map \eqn{D} that selects or combines coefficients, a scalar function
#' \eqn{\rho} applied to what the map returns, and hyperparameters \eqn{\theta}
#' that scale or shape it. The class holds no mathematics of its own. It
#' records the pieces that every branch needs, so that a penalty's size, its
#' hyperparameter names and their bounds can be read without knowing which
#' branch it is.
#'
#' @details
#' # The properties
#'
#' With \eqn{q} coefficients and a map of \eqn{m} rows, `beta` is a vector of
#' length \eqn{q}, \eqn{D\beta} has length \eqn{m}, and every derivative in the
#' coefficients is returned with length \eqn{q} or shape \eqn{q \times q}. A
#' `NULL` map is the identity, and then \eqn{m = q} and no arithmetic is done.
#'
#' Each hyperparameter carries a \pkg{linkfunctions7} link in `link_params`,
#' which maps its own open interval onto the whole real line. A penalty can
#' therefore be optimized on the unconstrained scale: a caller works in
#' \eqn{\eta = g(\theta)}, without bound constraints, and
#' [penalty_grad_theta()], [penalty_hess_theta()] and [penalty_cross()]
#' return either scale through their `scale` argument.
#'
#' `params_smooth` records in which hyperparameters the value is
#' differentiable. It is `TRUE` for every hyperparameter of the quadratic,
#' structured, additive, SCAD and MCP branches. A separable penalty copies the
#' entries of its parent distribution, so a penalty built by
#' [distrib_penalty()] on a Laplace distribution with a free location has
#' `FALSE` for that location.
#'
#' # Writing a branch
#'
#' This class is constructed directly only to write a new branch. A subclass
#' supplies [penalty_value()], [penalty_gradient()], [penalty_hessian()],
#' [penalty_grad_theta()], [penalty_hess_theta()], [penalty_cross()],
#' [penalty_kinks()] and [is_proper()], which have no method on the base
#' class, and [check_penalty()] then compares the derivatives with numerical
#' derivatives of the value. The methods registered on the base class, listed
#' below, supply defaults for everything else.
#'
#' The shipped branches are [quadratic_penalty()] (with [ridge_penalty()]),
#' [distrib_penalty()] (with [lasso_penalty()] and [elasticnet_penalty()]),
#' [scad_penalty()], [mcp_penalty()], [additive_penalty()] and
#' [structured_penalty()], the last built on a \pkg{parameters7} matrix
#' parameter.
#'
#' @section Methods:
#' The methods registered on `penalty` itself are inherited by every branch
#' that does not register its own.
#'
#' Predicates and readable quantities:
#' * [`is_quadratic()`][is_quadratic.penalty] returns `FALSE`.
#' * [`penalty_readable()`][penalty_readable.penalty] returns `NULL`.
#' * [`beta_quadratic()`][beta_quadratic.penalty] returns `FALSE`.
#'
#' Quantities of a quadratic penalty, which signal an error:
#' * [`penalty_matrix()`][penalty_matrix.penalty],
#'   [`penalty_rank()`][penalty_rank.penalty],
#'   [`penalty_null_basis()`][penalty_null_basis.penalty] and
#'   [`penalty_logpdet()`][penalty_logpdet.penalty].
#'
#' Derivatives for a marginal criterion, which signal an error:
#' * [`penalty_dhessian()`][penalty_dhessian.penalty],
#'   [`penalty_d2hessian()`][penalty_d2hessian.penalty] and
#'   [`penalty_dcross()`][penalty_dcross.penalty];
#' * [`penalty_dhessian_beta()`][penalty_dhessian_beta.penalty];
#' * [`penalty_d2hessian_beta()`][penalty_d2hessian_beta.penalty] and
#'   [`penalty_dhessian_beta_theta()`][penalty_dhessian_beta_theta.penalty].
#'
#' Proximal operator:
#' * [`penalty_prox()`][penalty_prox.penalty] signals an error.
#' * [`penalty_prox_spec()`][penalty_prox_spec.penalty] returns `NULL`.
#'
#' Draws and printing:
#' * [`penalty_draw()`][penalty_draw.penalty] draws from a quadratic penalty
#'   read as a Gaussian prior and returns `NA` otherwise.
#' * [`print()`][print.penalty] writes one line describing the penalty.
#'
#' @param penalty_name A single string naming the penalty, used by `print()`
#'   and by any code that reports which penalty a block carries.
#' @param map The matrix \eqn{D}, of \eqn{m} rows and \eqn{q} columns, or
#'   `NULL` for the identity. A \pkg{Matrix} object is kept in its own
#'   storage, so a diagonal map, which is what standardization produces, is
#'   recognized as diagonal by its class.
#' @param n_coef The number of coefficients \eqn{q}. A single whole number.
#' @param params The hyperparameter names, in the order every derivative list
#'   is keyed by. `character(0)` for a penalty with none.
#' @param params_bounds A named list, one entry per hyperparameter, each a
#'   numeric pair giving an **open** interval. A value at either endpoint is
#'   rejected, so `(0, Inf)` excludes zero.
#' @param link_params A named list, one \pkg{linkfunctions7} link per
#'   hyperparameter, carrying that hyperparameter's own interval onto the whole
#'   real line.
#' @param params_smooth A logical vector, one entry per hyperparameter, `TRUE`
#'   where the value is differentiable in it.
#'
#' @return An S7 object of class `penalty` carrying the seven properties
#'   above. The class is abstract: an object of exactly this class can be
#'   printed, it receives the defaults listed under Methods, and
#'   [penalty_value()] and the other generics without a base method signal an
#'   error for it.
#'
#' @seealso [quadratic_penalty()], [distrib_penalty()], [scad_penalty()],
#'   [mcp_penalty()], [additive_penalty()] and [structured_penalty()] for the
#'   branches;
#'   [penalty_value()] and [penalty_gradient()] for what a branch supplies;
#'   [check_penalty()] to verify one.
#'
#' @examples
#' # Every branch inherits from this class.
#' pen <- quadratic_penalty(crossprod(diff(diag(4))), map = NULL)
#' S7::S7_inherits(pen, penalty)
#'
#' # The properties that can be read without knowing the branch.
#' pen@penalty_name
#' pen@n_coef
#' pen@params
#' pen@params_bounds
#'
#' # A first-difference map over three coefficients: three coefficients in,
#' # two rows out.
#' D <- diff(diag(3))
#' mapped <- quadratic_penalty(diag(2), map = D)
#' mapped
#'
#' @export
penalty <- S7::new_class(
  name = "penalty",
  properties = list(
    penalty_name = S7::class_character,
    map = S7::class_any,
    n_coef = S7::class_numeric,
    params = S7::class_character,
    params_bounds = S7::class_list,
    link_params = S7::class_list,
    params_smooth = S7::class_logical
  )
)

#' Align and Validate the Hyperparameters
#'
#' @description
#' Puts a hyperparameter argument into the single shape that every branch
#' reads: a list reordered to `pen@params`, with the names removed from the
#' values, and checked against `pen@params_bounds` as open intervals. Returns
#' the aligned list. For a penalty with no hyperparameters it returns an empty
#' list without examining `theta`.
#'
#' @details
#' A named numeric vector is accepted in place of a list and converted with
#' `as.list()`, because the SCAD and MCP methods read their hyperparameters
#' with `$`, which does not apply to an atomic vector. The conversion is made
#' here, in the function that every generic calls, so every branch receives a
#' list.
#'
#' The bounds are **open**, so a hyperparameter at an endpoint is rejected
#' rather than clamped: `alpha = 1` on an elastic net whose bound is
#' \eqn{(0, 1)} signals an error. Parameters in \pkg{distributions7} are
#' validated in the same way, and the open bounds keep the value of every
#' link finite.
#'
#' @param pen A [penalty()] object.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same. Extra entries are dropped; a missing one is an
#'   error naming which. Each value must be a single number.
#'
#' @return A list of the same length and order as `pen@params`, each element
#'   an unnamed single number. `list()` when the penalty has no
#'   hyperparameters.
#'
#' @section Errors:
#' `Missing parameter(s) in 'theta': ... Expected: ...` when a name is absent
#' or `theta` is unnamed,
#' `Parameter 'p' must be a single number.` when a value is `NULL`, empty or
#' longer than one, and
#' `Parameter 'p' must lie in the open interval (a, b).` when a value is
#' non-finite or outside its bounds.
#'
#' @keywords internal
align_ptheta <- function(pen, theta) {
  params <- pen@params
  if (length(params) == 0L) return(list())
  # A named numeric vector carries what the list carries, and the branches
  # split on how they read it: `[[` accepts both, `$` accepts only the list.
  # A caller passing a vector therefore reached the quadratic and separable
  # branches and failed inside scad and mcp, three frames down and naming
  # neither the argument nor the penalty. The shape is settled here, at the
  # one point every generic passes through.
  if (!is.list(theta)) theta <- as.list(theta)
  if (is.null(names(theta)) || !all(params %in% names(theta))) {
    stop(sprintf("Missing parameter(s) in 'theta': %s. Expected: %s.",
                 paste(setdiff(params, names(theta)), collapse = ", "),
                 paste(params, collapse = ", ")), call. = FALSE)
  }
  theta <- theta[params]
  for (p in params) {
    b <- pen@params_bounds[[p]]
    v <- theta[[p]]
    if (length(v) != 1L) {
      stop(sprintf("Parameter '%s' must be a single number.", p),
           call. = FALSE)
    }
    if (any(!is.finite(v)) || any(v <= b[1]) || any(v >= b[2])) {
      stop(sprintf("Parameter '%s' must lie in the open interval (%g, %g).",
                   p, b[1], b[2]), call. = FALSE)
    }
    theta[[p]] <- unname(v)
  }
  theta
}

#' A Map, in Whatever Storage It Arrived In
#'
#' @description
#' Returns the map unchanged when it is already a \pkg{Matrix} object, and
#' `as.matrix()` of it otherwise. Called once, by each constructor that takes
#' a map, so that a map given as a data frame or a vector becomes a matrix
#' while a sparse or diagonal one keeps its own storage.
#'
#' @details
#' A dense copy of a diagonal map would hold \eqn{q^2} numbers instead of
#' \eqn{q}. A diagonal map is what standardization produces, a rescaling of
#' each coordinate under which a separable penalty stays separable and its
#' proximal operator stays in closed form. The product and the crossproduct,
#' which are the only operations applied to the map, are defined for both
#' kinds of storage.
#'
#' @param map A matrix, a \pkg{Matrix}, or anything `as.matrix()` accepts.
#'
#' @return The same object when it is a \pkg{Matrix}, and a base matrix
#'   otherwise.
#'
#' @keywords internal
as_map <- function(map) {
  if (isS4(map) && methods::is(map, "Matrix")) return(map)
  as.matrix(map)
}

#' Apply the Linear Map and Its Transpose
#'
#' @description
#' `map_apply()` computes \eqn{t = D\beta}, carrying a coefficient vector to
#' the argument at which \eqn{\rho} is evaluated. `map_back()` computes
#' \eqn{D'g},
#' carrying a gradient in \eqn{t} back to a gradient in \eqn{\beta}. A `NULL`
#' map is the identity and both return their argument untouched.
#'
#' @param pen A [penalty()] object, whose `map` is \eqn{D} with \eqn{m} rows
#'   and \eqn{q} columns, or `NULL`.
#' @param beta A numeric vector of length \eqn{q}. `map_apply()` only.
#' @param g A numeric vector of length \eqn{m}. `map_back()` only.
#'
#' @return `map_apply()` a numeric vector of length \eqn{m}; `map_back()` a
#'   numeric vector of length \eqn{q}. Both are plain numeric vectors even
#'   when the map is a \pkg{Matrix} object, rather than one-column
#'   \pkg{Matrix} objects.
#'
#' @keywords internal
map_apply <- function(pen, beta) {
  if (is.null(pen@map)) beta else as.numeric(pen@map %*% beta)
}

#' @rdname map_apply
#' @keywords internal
map_back <- function(pen, g) {
  if (is.null(pen@map)) g else as.numeric(crossprod(pen@map, g))
}

#' Carry a Middle Matrix Through the Map
#'
#' @description
#' `map_quad()` computes \eqn{D' \mathrm{diag}(h) D} without forming the
#' diagonal matrix, which is the Hessian of a separable penalty carried back to
#' the coefficients. `map_quad_full()` computes \eqn{D'MD} for a parent read
#' blockwise, whose middle matrix is block diagonal. With a `NULL` map the
#' first returns `diag(h)` and the second returns \eqn{M}.
#'
#' @details
#' Both return a base matrix even when the map is a \pkg{Matrix} object. The
#' crossproduct of a \pkg{Matrix} map would otherwise be a \pkg{Matrix}
#' object, whereas the result for a `NULL` map is a base matrix and
#' [map_back()] returns a plain vector. The conversion is made here so that
#' the caller receives a base matrix from every branch and can write it into
#' a block of a base information matrix.
#'
#' @param pen A [penalty()] object, whose `map` is \eqn{D} with \eqn{m} rows
#'   and \eqn{q} columns, or `NULL`.
#' @param h A numeric vector of length \eqn{m}, the diagonal entries.
#'   `map_quad()` only.
#' @param m A symmetric \eqn{m \times m} matrix. `map_quad_full()` only.
#'
#' @return A \eqn{q \times q} symmetric base matrix, from both.
#'
#' @keywords internal
map_quad <- function(pen, h) {
  if (is.null(pen@map)) return(diag(h, length(h)))
  # a Matrix map carries its class through the crossproduct, and the result
  # would then be the one thing in the contract that is not a base matrix:
  # the identity-map branch above is already dense at any width, map_back()
  # coerces its vector for the same reason, and a consumer writing this into
  # a block of its own information fails on the class rather than on the
  # arithmetic. Densified here, where the contract is stated.
  as.matrix(crossprod(pen@map, pen@map * h))
}

#' @rdname map_quad
#' @keywords internal
map_quad_full <- function(pen, m) {
  if (is.null(pen@map)) return(as.matrix(m))
  as.matrix(crossprod(pen@map, m %*% pen@map))
}

#' The Hyperparameter Pair Names
#'
#' @description
#' The keys of a penalty's second hyperparameter derivatives, and the pairs
#' they stand for: the \eqn{p} diagonals first, in `params` order, then the
#' \eqn{p(p-1)/2} upper off-diagonal pairs, each joined by an underscore. For
#' `c("lambda", "alpha")` the keys are `lambda_lambda`, `alpha_alpha`,
#' `lambda_alpha`. This is the order [penalty_hess_theta()] returns and the
#' order [ptheta_to_link()] reads.
#'
#' @details
#' The diagonals come first, rather than in lexicographic order, so that the
#' first \eqn{p} entries are the diagonal of the Hessian. The same convention
#' names the Hessian components of \pkg{distributions7}.
#'
#' `character(0)` gives an empty **named** list. A penalty with no free
#' hyperparameters has no pair in which to differentiate, and the named empty
#' list is the shape that [penalty_grad_theta()] and [penalty_cross()] return
#' in that case. The function treats this case separately, because
#' `paste0(character(0), "_", character(0))` returns the single string `"_"`
#' and [stats::setNames()] would then signal an error on the empty list of
#' pairs.
#'
#' @param params A character vector of hyperparameter names, in the order in
#'   which the penalty holds them. May be empty.
#'
#' @return A named list of length \eqn{p(p+1)/2}, empty when `params` is. Each
#'   element is a character pair naming the two hyperparameters of the
#'   derivative, and each name is those two joined by an underscore.
#'
#' @keywords internal
ptheta_pairs <- function(params) {
  p <- length(params)
  # paste0 recycles a zero-length argument against the length-one literal, so
  # without this the names are the single string "_" while the list of pairs is
  # empty, and setNames() raises. A penalty with no free hyperparameters has no
  # pairs to differentiate in; align_ptheta() answers the same way.
  if (p == 0L) return(stats::setNames(list(), character(0)))
  nm <- paste0(params, "_", params)
  prs <- lapply(params, function(x) c(x, x))
  if (p > 1L) {
    for (i in seq_len(p - 1L)) {
      for (j in seq.int(i + 1L, p)) {
        nm <- c(nm, paste0(params[i], "_", params[j]))
        prs <- c(prs, list(c(params[i], params[j])))
      }
    }
  }
  stats::setNames(prs, nm)
}

#' Carry Hyperparameter Derivatives Onto the Unconstrained Scale
#'
#' @description
#' Applies the chain rule that turns a derivative in \eqn{\theta} into one in
#' \eqn{\eta = g(\theta)}, at first and second order, using the diagonal
#' Jacobian the links supply. Handles one of the three derivative kinds per
#' call, whichever of `g`, `H` and `cross` is given.
#'
#' @details
#' The links are scalar and one per hyperparameter, so the Jacobian is diagonal
#' and the chain rule needs no partition sums. Writing \eqn{h = g^{-1}} and
#' \eqn{\eta_i = g_i(\theta_i)},
#'
#' \deqn{\frac{\partial\rho}{\partial\eta_i}
#'   = \frac{\partial\rho}{\partial\theta_i}\, h_i'(\eta_i), \qquad
#'   \frac{\partial^2\rho}{\partial\eta_i \partial\eta_j}
#'   = \frac{\partial^2\rho}{\partial\theta_i \partial\theta_j}\,
#'     h_i'(\eta_i)\, h_j'(\eta_j)
#'   + \delta_{ij}\, \frac{\partial\rho}{\partial\theta_i}\, h_i''(\eta_i).}
#'
#' The second-derivative term appears on the diagonal alone, which is why the
#' Hessian branch requires the gradient as well. The mixed block
#' \eqn{\partial^2\rho / \partial\beta\,\partial\theta_i} is first order in
#' \eqn{\theta} and picks up one factor of \eqn{h_i'}, because a
#' reparametrization of the hyperparameters does not change the coefficient
#' direction.
#'
#' This is the chain rule that \pkg{distributions7} applies to the link scale
#' of a distribution's parameters, restricted to the first two orders.
#'
#' @param pen A [penalty()] object.
#' @param theta The aligned hyperparameters, as [align_ptheta()] returns them.
#' @param g The parameter-scale gradient list, keyed by `pen@params`. Required
#'   for the gradient and for the Hessian; `NULL` otherwise.
#' @param H The parameter-scale Hessian list, keyed as [ptheta_pairs()]
#'   names the pairs. Supplied with `g` for the second order.
#' @param cross The parameter-scale mixed list, keyed by `pen@params`, each
#'   element a vector of length `pen@n_coef`. Supplied alone.
#'
#' @return Whichever kind was supplied, on the unconstrained scale, with the
#'   same names and shapes as the input. `cross` takes precedence over `H`,
#'   and `H` over `g`, when more than one is given.
#'
#' @keywords internal
ptheta_to_link <- function(pen, theta, g = NULL, H = NULL, cross = NULL) {
  params <- pen@params
  h1 <- h2 <- stats::setNames(numeric(length(params)), params)
  for (p in params) {
    lk <- pen@link_params[[p]]
    eta <- linkfunctions7::linkfun(lk, theta[[p]])
    h1[p] <- linkfunctions7::dlinkinv(lk, eta)
    h2[p] <- linkfunctions7::d2linkinv(lk, eta)
  }
  if (!is.null(cross)) {
    return(stats::setNames(lapply(params, function(p) cross[[p]] * h1[[p]]),
                           params))
  }
  if (is.null(H)) {
    return(stats::setNames(lapply(params, function(p) g[[p]] * h1[[p]]),
                           params))
  }
  prs <- ptheta_pairs(params)
  stats::setNames(lapply(names(prs), function(nm) {
    pr <- prs[[nm]]
    out <- H[[nm]] * h1[[pr[1]]] * h1[[pr[2]]]
    if (pr[1] == pr[2]) out <- out + g[[pr[1]]] * h2[[pr[1]]]
    out
  }), names(prs))
}

#' @title Print a Penalty
#' @name print.penalty
#'
#' @description
#' Writes one line giving the name of the penalty, the number of coefficients,
#' the number of rows of its map, and its hyperparameters. The row count is
#' `nrow(map)`, or `n_coef` when the map is `NULL`, so an unmapped penalty
#' shows the same number twice.
#'
#' @param x A [penalty()] object of any branch.
#' @param ... Unused, accepted for consistency with [print()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' quadratic_penalty(diag(3))
#'
#' # A first-difference map: three coefficients, two differences.
#' quadratic_penalty(diag(2), map = diff(diag(3)))
#'
#' # A penalty with no hyperparameters prints "(none)".
#' distrib_penalty(
#'   distributions7::fixed(distributions7::gaussian1_distrib(),
#'                         mu = 0, sigma = 1),
#'   n_coef = 3)
#'
#' @keywords internal
S7::method(print, penalty) <- function(x, ...) {
  m <- if (is.null(x@map)) x@n_coef else nrow(x@map)
  cat(sprintf("%s penalty on %d coefficient(s) through %d row(s); theta: %s\n",
              x@penalty_name, as.integer(x@n_coef), as.integer(m),
              if (length(x@params)) paste(x@params, collapse = ", ")
              else "(none)"))
  invisible(x)
}
