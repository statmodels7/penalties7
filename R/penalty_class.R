#' @title S7 Base Class for Penalties
#'
#' @description
#' The abstract parent of every penalty. A penalty is
#' \eqn{\rho(D\beta;\theta)}: a linear map \eqn{D}, a scalar function
#' \eqn{\rho} and hyperparameters \eqn{\theta}. Each hyperparameter travels
#' with a \pkg{linkfunctions7} link, so a consumer can optimize on the
#' unconstrained scale.
#'
#' @param penalty_name A string naming the penalty.
#' @param map The matrix \eqn{D}, or \code{NULL} for the identity.
#' @param n_coef The number of coefficients \eqn{q}.
#' @param params Hyperparameter names, in order.
#' @param params_bounds A named list of open intervals.
#' @param link_params A named list of \pkg{linkfunctions7} links.
#' @param params_smooth Logical vector; which hyperparameters are
#'   differentiable.
#'
#' @return An object inheriting from class \code{penalty}.
#'
#' @seealso \code{\link{quadratic_penalty}}, \code{\link{distrib_penalty}},
#'   \code{\link{scad_penalty}}
#'
#' @examples
#' S7::S7_inherits(quadratic_penalty(diag(3)), penalty)
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
#' Reorders \code{theta} by name, strips stray names off the values and
#' validates against \code{params_bounds} treated as open intervals -- the
#' \pkg{distributions7} contract, restated here for hyperparameters. A named
#' numeric vector is accepted in place of the list and converted to one, so
#' that every branch reads the same shape.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same.
#'
#' @return The aligned list.
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
    if (any(!is.finite(v)) || any(v <= b[1]) || any(v >= b[2])) {
      stop(sprintf("Parameter '%s' must lie in the open interval (%g, %g).",
                   p, b[1], b[2]), call. = FALSE)
    }
    theta[[p]] <- unname(v)
  }
  theta
}

#' A Map, in Whatever Form It Keeps
#'
#' @description
#' The map as the caller gave it, densified only where it is not already a
#' matrix of some kind.
#'
#' @details
#' A \pkg{Matrix} object is kept as it is. Densifying a diagonal map would
#' cost \eqn{q^2} numbers where it holds \eqn{q}, and a diagonal map is
#' exactly what standardization is: a rescaling of each coordinate, under
#' which a separable penalty stays separable and its proximal operator stays
#' closed. Every arithmetic the map takes part in -- the product, the
#' crossproduct -- works for both kinds.
#'
#' @param map A matrix, a \pkg{Matrix}, or anything coercible to one.
#'
#' @return The map.
#'
#' @keywords internal
as_map <- function(map) {
  if (isS4(map) && methods::is(map, "Matrix")) return(map)
  as.matrix(map)
}

#' Apply the Linear Map and Its Transpose
#'
#' @description
#' \code{map_apply} computes \eqn{t = D\beta} and \code{map_back} computes
#' \eqn{D'g}; a \code{NULL} map is the identity and pays nothing.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param beta A numeric vector of coefficients.
#' @param g A numeric vector of length \code{nrow(D)}.
#'
#' @return A numeric vector.
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

#' Carry a Diagonal Middle Matrix Through the Map
#'
#' @description
#' \eqn{D' \mathrm{diag}(h) D} for the separable Hessians, without forming
#' the diagonal matrix.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param h A numeric vector of diagonal entries.
#'
#' @return A \code{q x q} symmetric matrix.
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

#' The Hyperparameter Pair Names
#'
#' @description
#' The component names of a penalty's second theta derivatives: diagonals
#' first, then the upper off-diagonal pairs, joined by an underscore.
#'
#' @param params The hyperparameter names.
#'
#' @return A character vector.
#'
#' @keywords internal
ptheta_pairs <- function(params) {
  p <- length(params)
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

#' Carry Theta Derivatives Onto the Link Scale
#'
#' @description
#' The order 1-2 chain rule with the diagonal Jacobian of the links: the
#' same interception \pkg{distributions7} applies, restricted to the two
#' orders a penalty consumer needs.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param theta The aligned hyperparameters.
#' @param g The parameter-scale gradient list, or \code{NULL}.
#' @param H The parameter-scale Hessian list, or \code{NULL}.
#' @param cross The parameter-scale mixed list, or \code{NULL}.
#'
#' @return Whichever of the three was supplied, transformed.
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
#' @description One line: the name, the sizes, the hyperparameters.
#' @param x A \code{\link{penalty}} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @keywords internal
S7::method(print, penalty) <- function(x, ...) {
  m <- if (is.null(x@map)) x@n_coef else nrow(x@map)
  cat(sprintf("%s penalty on %d coefficient(s) through %d row(s); theta: %s\n",
              x@penalty_name, as.integer(x@n_coef), as.integer(m),
              if (length(x@params)) paste(x@params, collapse = ", ")
              else "(none)"))
  invisible(x)
}
