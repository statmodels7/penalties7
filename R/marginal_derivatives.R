#' @include generics.R quadratic_penalty.R additive_penalty.R structured_penalty.R distrib_penalty.R scad_mcp.R
NULL

# What a marginal criterion asks of a penalty, beyond the value, the gradient,
# the Hessian and the mixed block.
#
# Estimating the hyperparameters by a marginal likelihood needs the derivative
# of the LAPLACE approximation, whose determinant is |H + S| with S the
# penalty's Hessian in the coefficients. Differentiating it in the
# hyperparameters therefore needs dS/dtheta, and differentiating a second time
# needs d2S/dtheta2 and d3rho/dbeta dtheta2. None of those is expressible from
# the value, the gradient and the Hessian: they are new quantities, and they
# belong here rather than in whatever consumes them, so that a penalty written
# later answers them by writing three methods and works in a marginal
# criterion with no change anywhere else.
#
# The base class rejects. A penalty that cannot answer says so, rather than
# being read through a test of its behaviour -- a consumer that measured
# whether a Hessian happened to be linear in the hyperparameters would be
# guessing at a property the penalty knows.

#' The Derivative of the Coefficient Hessian in the Hyperparameters
#'
#' @description
#' \eqn{\partial S/\partial\theta_m = \partial^3\rho/\partial\beta^2\,
#' \partial\theta_m}, one matrix per hyperparameter.
#'
#' @details
#' \code{penalty_hessian} says how curved the penalty is in the coefficients;
#' this says how that curvature moves with each hyperparameter. It is the piece
#' a marginal criterion needs to differentiate \eqn{\log|H+S|}.
#'
#' Every branch answers in closed form. For \code{\link{quadratic_penalty}} the
#' Hessian is \eqn{\lambda D'PD} and the derivative is \eqn{D'PD}; for
#' \code{\link{additive_penalty}} it is the component \eqn{P_k} of each
#' smoothing parameter; for \code{\link{structured_penalty}} it is the matrix
#' parameter's own \code{param_d1}; for \code{\link{distrib_penalty}} it is
#' \eqn{-D'\mathrm{diag}(\partial^3\ell/\partial y^2\partial\theta_m)D}, which
#' \pkg{distributions7} supplies as \code{distrib_cross2_y}. A penalty defined
#' by its derivative rather than by a density -- \code{\link{scad_penalty}},
#' \code{\link{mcp_penalty}} -- has no such quantity where its kinks are and
#' rejects.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameters.
#' @param ... Passed to methods.
#'
#' @return A named list of matrices, one per hyperparameter.
#'
#' @seealso \code{\link{penalty_hessian}}, \code{\link{penalty_d2hessian}},
#'   \code{\link{penalty_dcross}}
#'
#' @examples
#' pen <- quadratic_penalty(diag(3))
#' penalty_dhessian(pen, c(1, 2, 3), list(lambda = 2))
#'
#' @export
penalty_dhessian <- S7::new_generic("penalty_dhessian", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    S7::S7_dispatch()
  })


#' The Second Derivative of the Coefficient Hessian in the Hyperparameters
#'
#' @description
#' \eqn{\partial^2 S/\partial\theta_m\partial\theta_l =
#' \partial^4\rho/\partial\beta^2\partial\theta_m\partial\theta_l}, one matrix
#' per unordered pair.
#'
#' @details
#' Zero for every penalty whose Hessian is linear in its hyperparameters, which
#' is the quadratic and additive branches; the structured branch reads the
#' matrix parameter's \code{param_d2}, and the separable branch the second
#' \eqn{\theta}-derivative of the parent's response curvature.
#'
#' The keys are those of \code{\link{penalty_hess_theta}}, so a consumer
#' looking a pair up need not know which order it was written in.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameters.
#' @param ... Passed to methods.
#'
#' @return A named list of matrices, keyed by hyperparameter pair.
#'
#' @seealso \code{\link{penalty_dhessian}}
#'
#' @examples
#' pen <- quadratic_penalty(diag(3))
#' penalty_d2hessian(pen, c(1, 2, 3), list(lambda = 2))
#'
#' @export
penalty_d2hessian <- S7::new_generic("penalty_d2hessian", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    S7::S7_dispatch()
  })


#' The Derivative of the Mixed Block in the Hyperparameters
#'
#' @description
#' \eqn{\partial^3\rho/\partial\beta\,\partial\theta_m\partial\theta_l}, one
#' vector per unordered pair.
#'
#' @details
#' \code{\link{penalty_cross}} is how the coefficient gradient moves with one
#' hyperparameter; this is how that movement itself moves with a second. It is
#' what a marginal criterion needs to differentiate the mode's own derivative,
#' and it is exactly zero wherever the penalty is quadratic in the
#' coefficients with a Hessian linear in the hyperparameters.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameters.
#' @param ... Passed to methods.
#'
#' @return A named list of numeric vectors, keyed by hyperparameter pair.
#'
#' @seealso \code{\link{penalty_cross}}, \code{\link{penalty_dhessian}}
#'
#' @examples
#' pen <- quadratic_penalty(diag(3))
#' penalty_dcross(pen, c(1, 2, 3), list(lambda = 2))
#'
#' @export
penalty_dcross <- S7::new_generic("penalty_dcross", "pen",
  function(pen, beta, theta, ...) {
    theta <- align_ptheta(pen, theta)
    beta <- as.numeric(beta)
    S7::S7_dispatch()
  })


#' Is a Penalty Quadratic in the Coefficients?
#'
#' @description
#' \code{TRUE} when \eqn{\partial^3\rho/\partial\beta^3} is exactly zero, which
#' spares a consumer that third derivative altogether.
#'
#' @details
#' Distinct from \code{\link{is_quadratic}}, which is about the whole
#' construction of \code{\link{quadratic_penalty}}: a structured penalty is
#' quadratic in the coefficients and is not a quadratic penalty, and a
#' \code{\link{distrib_penalty}} over a gaussian parent is quadratic in the
#' coefficients while its Hessian is not linear in its hyperparameters. The
#' three properties are independent and each is asked of the penalty.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param theta A named list of hyperparameters.
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{penalty_dhessian}}
#'
#' @examples
#' beta_quadratic(quadratic_penalty(diag(3)), list(lambda = 1))
#'
#' @export
beta_quadratic <- S7::new_generic("beta_quadratic", "pen",
  function(pen, theta, ...) S7::S7_dispatch())


# --- the base class rejects -------------------------------------------------

#' @title What the Base Class Answers
#' @name penalty_dhessian.penalty
#' @description
#' Each of these rejects, naming what is missing. A marginal criterion is a
#' Laplace approximation and asks for derivatives beyond the second; a penalty
#' that does not supply them cannot be estimated by one, and reporting that is
#' better than a criterion assembled from a quantity nobody wrote.
#' @param pen A \code{\link{penalty}} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameters.
#' @param ... Unused.
#' @return Signals an error.
#' @keywords internal
S7::method(penalty_dhessian, penalty) <- function(pen, beta, theta, ...) {
  stop(sprintf(paste0("'%s' does not supply penalty_dhessian(), so its\n",
                      "  hyperparameters cannot be estimated by a marginal",
                      " criterion."),
               pen@penalty_name), call. = FALSE)
}

#' @rdname penalty_dhessian.penalty
#' @name penalty_d2hessian.penalty
#' @keywords internal
S7::method(penalty_d2hessian, penalty) <- function(pen, beta, theta, ...) {
  stop(sprintf("'%s' does not supply penalty_d2hessian().",
               pen@penalty_name), call. = FALSE)
}

#' @rdname penalty_dhessian.penalty
#' @name penalty_dcross.penalty
#' @keywords internal
S7::method(penalty_dcross, penalty) <- function(pen, beta, theta, ...) {
  stop(sprintf("'%s' does not supply penalty_dcross().",
               pen@penalty_name), call. = FALSE)
}

#' @rdname penalty_dhessian.penalty
#' @name beta_quadratic.penalty
#' @keywords internal
S7::method(beta_quadratic, penalty) <- function(pen, theta, ...) FALSE


# --- quadratic --------------------------------------------------------------

#' @title Marginal Derivatives of a Quadratic Penalty
#' @name penalty_dhessian.QuadraticPenalty
#' @description
#' \eqn{S = \lambda D'PD} is linear in \eqn{\lambda} and free of the
#' coefficients, so \eqn{\partial S/\partial\lambda = D'PD} and every higher
#' derivative is zero.
#' @param pen A \code{QuadraticPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameters.
#' @param ... Unused.
#' @return A named list of matrices, of vectors, or a logical.
#' @keywords internal
S7::method(penalty_dhessian, QuadraticPenalty) <- function(pen, beta, theta,
                                                           ...) {
  # `unclass()` strips the attributes a base matrix may carry and does
  # nothing to an S4 one, which is what a blocked penalty stores: the
  # derivative is the matrix itself, in whatever storage the penalty keeps
  # it, and a caller that needs a base matrix coerces where the two meet.
  list(lambda = if (isS4(pen@DPD)) pen@DPD else unclass(pen@DPD))
}

#' @rdname penalty_dhessian.QuadraticPenalty
#' @name penalty_d2hessian.QuadraticPenalty
#' @keywords internal
S7::method(penalty_d2hessian, QuadraticPenalty) <- function(pen, beta, theta,
                                                            ...) {
  zero_pairs(pen, matrix(0, as.integer(pen@n_coef), as.integer(pen@n_coef)))
}

#' @rdname penalty_dhessian.QuadraticPenalty
#' @name penalty_dcross.QuadraticPenalty
#' @keywords internal
S7::method(penalty_dcross, QuadraticPenalty) <- function(pen, beta, theta,
                                                         ...) {
  zero_pairs(pen, numeric(as.integer(pen@n_coef)))
}

#' @rdname penalty_dhessian.QuadraticPenalty
#' @name beta_quadratic.QuadraticPenalty
#' @keywords internal
S7::method(beta_quadratic, QuadraticPenalty) <- function(pen, theta, ...) TRUE


# --- additive ---------------------------------------------------------------

#' @title Marginal Derivatives of an Additive Penalty
#' @name penalty_dhessian.AdditivePenalty
#' @description
#' \eqn{S = \sum_k \lambda_k P_k} is linear in the smoothing parameters, so
#' \eqn{\partial S/\partial\lambda_k = P_k} and every higher derivative is
#' zero.
#' @param pen An \code{AdditivePenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameters.
#' @param ... Unused.
#' @return A named list of matrices, of vectors, or a logical.
#' @keywords internal
S7::method(penalty_dhessian, AdditivePenalty) <- function(pen, beta, theta,
                                                          ...) {
  stats::setNames(lapply(pen@mats, function(P) unclass(P)), pen@params)
}

#' @rdname penalty_dhessian.AdditivePenalty
#' @name penalty_d2hessian.AdditivePenalty
#' @keywords internal
S7::method(penalty_d2hessian, AdditivePenalty) <- function(pen, beta, theta,
                                                           ...) {
  zero_pairs(pen, matrix(0, as.integer(pen@n_coef), as.integer(pen@n_coef)))
}

#' @rdname penalty_dhessian.AdditivePenalty
#' @name penalty_dcross.AdditivePenalty
#' @keywords internal
S7::method(penalty_dcross, AdditivePenalty) <- function(pen, beta, theta,
                                                        ...) {
  zero_pairs(pen, numeric(as.integer(pen@n_coef)))
}

#' @rdname penalty_dhessian.AdditivePenalty
#' @name beta_quadratic.AdditivePenalty
#' @keywords internal
S7::method(beta_quadratic, AdditivePenalty) <- function(pen, theta, ...) TRUE


# --- structured -------------------------------------------------------------

#' @title Marginal Derivatives of a Structured Penalty
#' @name penalty_dhessian.StructuredPenalty
#' @description
#' \eqn{S = \Omega(\theta)} is the matrix parameter itself, so its derivatives
#' in the hyperparameters are the structure's own \code{param_d1} and
#' \code{param_d2}, which \pkg{parameters7} supplies exactly. It is quadratic
#' in the coefficients, so the mixed third derivative is
#' \eqn{\partial^2\Omega/\partial\theta_m\partial\theta_l\,\beta}.
#' @param pen A \code{StructuredPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameters.
#' @param ... Unused.
#' @return A named list of matrices, of vectors, or a logical.
#' @keywords internal
S7::method(penalty_dhessian, StructuredPenalty) <- function(pen, beta, theta,
                                                            ...) {
  d1 <- parameters7::param_d1(pen@structure, struct_eta(pen, theta))
  stats::setNames(lapply(seq_along(pen@params), function(m)
    unname(unclass(d1[[m]]))), pen@params)
}

#' @rdname penalty_dhessian.StructuredPenalty
#' @name penalty_d2hessian.StructuredPenalty
#' @keywords internal
S7::method(penalty_d2hessian, StructuredPenalty) <- function(pen, beta, theta,
                                                             ...) {
  d2 <- parameters7::param_d2(pen@structure, struct_eta(pen, theta))
  # param_d2 is a FLAT list, not a nested one, and its keys join the free
  # names with a separator of its own. The entry is located through
  # param_tuple_indices(), which enumerates the tuples in exactly the order
  # the names are in, so nothing here depends on the spelling of a key.
  idx <- parameters7::param_tuple_indices(pen@structure, 2L)
  prs <- ptheta_pairs(pen@params)
  stats::setNames(lapply(names(prs), function(nm) {
    ij <- sort(match(prs[[nm]], pen@params))
    hit <- which(vapply(idx, function(t) identical(sort(t), ij), logical(1)))
    if (!length(hit)) {
      stop(sprintf("No second-derivative component for '%s' and '%s'.",
                   prs[[nm]][1L], prs[[nm]][2L]), call. = FALSE)
    }
    unname(unclass(d2[[hit[1L]]]))
  }), names(prs))
}

#' @rdname penalty_dhessian.StructuredPenalty
#' @name penalty_dcross.StructuredPenalty
#' @keywords internal
S7::method(penalty_dcross, StructuredPenalty) <- function(pen, beta, theta,
                                                          ...) {
  d2 <- penalty_d2hessian(pen, beta, theta)
  lapply(d2, function(M) as.numeric(M %*% beta))
}

#' @rdname penalty_dhessian.StructuredPenalty
#' @name beta_quadratic.StructuredPenalty
#' @keywords internal
S7::method(beta_quadratic, StructuredPenalty) <- function(pen, theta, ...) TRUE


# --- separable --------------------------------------------------------------

#' @title Marginal Derivatives of a Separable Penalty
#' @name penalty_dhessian.DistribPenalty
#' @description
#' With \eqn{\rho = -\sum_j \log f((D\beta)_j;\theta)} the Hessian is
#' \eqn{-D'\mathrm{diag}(\ell^{(yy)})D}, so its \eqn{\theta}-derivatives are
#' the parent's \code{distrib_cross2_y} carried through the same map, and the
#' mixed third derivative is \code{distrib_cross_y} differentiated once more in
#' \eqn{\theta}. Both come from \pkg{distributions7} rather than being written
#' again here.
#' @details
#' The second derivatives read the parent's
#' \code{\link[distributions7]{distrib_grad_y_hess}} and
#' \code{\link[distributions7]{distrib_hess_y_hess}}, so nothing is
#' differentiated here either. A parent with closed forms for those -- the
#' gaussian, hence every ridge and every random effect -- makes this branch
#' exact; one without them inherits that package's documented fallback, which
#' is one central difference of its analytic first-order component.
#' @param pen A \code{DistribPenalty} object.
#' @param beta A numeric vector of coefficients.
#' @param theta A named list of hyperparameters.
#' @param ... Unused.
#' @return A named list of matrices, of vectors, or a logical.
#' @keywords internal
S7::method(penalty_dhessian, DistribPenalty) <- function(pen, beta, theta,
                                                         ...) {
  reject_kinked(pen, "penalty_dhessian")
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    c2 <- distributions7::distrib_cross2_y(pen@parent, a, theta)
    return(stats::setNames(lapply(pen@params, function(m)
      -map_quad_full(pen, dp_blockdiag(pen, c2[[m]], nrow(a)))), pen@params))
  }
  c2 <- distributions7::distrib_cross2_y(pen@parent, t, theta)
  stats::setNames(lapply(pen@params, function(m)
    -map_quad(pen, c2[[m]] + 0 * t)), pen@params)
}

#' @rdname penalty_dhessian.DistribPenalty
#' @name penalty_d2hessian.DistribPenalty
#' @keywords internal
S7::method(penalty_d2hessian, DistribPenalty) <- function(pen, beta, theta,
                                                          ...) {
  reject_kinked(pen, "penalty_d2hessian")
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    h <- distributions7::distrib_hess_y_hess(pen@parent, a, theta)
    return(carry_pairs(pen, h, function(v)
      -map_quad_full(pen, dp_blockdiag(pen, v, nrow(a)))))
  }
  h <- distributions7::distrib_hess_y_hess(pen@parent, t, theta)
  carry_pairs(pen, h, function(v) -map_quad(pen, v + 0 * t))
}

#' @rdname penalty_dhessian.DistribPenalty
#' @name penalty_dcross.DistribPenalty
#' @keywords internal
S7::method(penalty_dcross, DistribPenalty) <- function(pen, beta, theta, ...) {
  reject_kinked(pen, "penalty_dcross")
  t <- map_apply(pen, beta)
  if (pen@block > 1L) {
    a <- dp_arg(pen, t)
    g <- distributions7::distrib_grad_y_hess(pen@parent, a, theta)
    return(carry_pairs(pen, g, function(v)
      -map_back(pen, dp_flat(pen, v + 0 * a))))
  }
  g <- distributions7::distrib_grad_y_hess(pen@parent, t, theta)
  carry_pairs(pen, g, function(v) -map_back(pen, v + 0 * t))
}

#' @rdname penalty_dhessian.DistribPenalty
#' @name beta_quadratic.DistribPenalty
#' @keywords internal
S7::method(beta_quadratic, DistribPenalty) <- function(pen, theta, ...) {
  # the log-density is quadratic in the response exactly when its second
  # derivative there does not depend on it. The question is put that way, and
  # not to a third derivative, because the second is analytic for almost every
  # family while the third is often a difference, whose noise no threshold
  # separates from a true zero.
  t <- c(-1.73, -0.29, 0.61, 2.04)
  h <- tryCatch(distributions7::distrib_hess_y(pen@parent, t,
                                               align_ptheta(pen, theta)),
                error = function(e) NULL)
  !is.null(h) && length(h) > 1L &&
    isTRUE(all.equal(as.numeric(h), rep(as.numeric(h)[1L], length(h)),
                     tolerance = 1e-12))
}


#' Carry a Parent's Paired Components Into Coefficient Space
#'
#' @description
#' Re-keys a \pkg{distributions7} component keyed by parameter pair into this
#' package's own keys, and places each through the penalty's map.
#'
#' @details
#' The two enumerations of pairs are built the same way from the same names, so
#' a key from one is a key of the other; it is looked up by name in both
#' orders rather than by position, since a hyperparameter whose own name
#' contains the separator would not survive being taken apart.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param comp The parent's components, keyed by parameter pair.
#' @param carry A function placing one component into coefficient space.
#'
#' @return A named list keyed by hyperparameter pair.
#'
#' @keywords internal
carry_pairs <- function(pen, comp, carry) {
  prs <- ptheta_pairs(pen@params)
  keys <- names(comp)
  stats::setNames(lapply(names(prs), function(nm) {
    pr <- prs[[nm]]
    key <- paste(pr[1], pr[2], sep = "_")
    if (!key %in% keys) key <- paste(pr[2], pr[1], sep = "_")
    if (!key %in% keys) {
      stop(sprintf("The parent has no component for '%s' and '%s'.",
                   pr[1], pr[2]), call. = FALSE)
    }
    carry(comp[[key]])
  }), names(prs))
}


#' A Zero Entry for Every Hyperparameter Pair
#'
#' @description
#' The answer of a penalty whose Hessian is linear in its hyperparameters.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param z The zero object, a matrix or a vector.
#'
#' @return A named list keyed by hyperparameter pair.
#'
#' @keywords internal
zero_pairs <- function(pen, z) {
  nm <- names(ptheta_pairs(pen@params))
  stats::setNames(rep(list(z), length(nm)), nm)
}


#' Reject a Kinked Penalty
#'
#' @description
#' Signals that a penalty with a kink has no derivative of this order.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param what The generic's name.
#'
#' @return \code{NULL}, invisibly, or signals an error.
#'
#' @keywords internal
reject_kinked <- function(pen, what) {
  # the property is read directly rather than through penalty_kinks(), which
  # takes a theta this is asked before there is one to pass
  k <- if ("kinks" %in% S7::prop_names(pen)) pen@kinks else numeric(0)
  if (length(k) && any(is.finite(unlist(k)))) {
    stop(sprintf(paste0("'%s' has a kink, so %s() does not exist there and\n",
                        "  its hyperparameters cannot be estimated by a",
                        " marginal criterion."),
                 pen@penalty_name, what), call. = FALSE)
  }
  invisible(NULL)
}
