#' @include generics.R distrib_penalty.R quadratic_penalty.R prox.R
NULL

#' @title Drawing From a Penalty Read as a Prior
#'
#' @description
#' A vector drawn from the distribution of which the penalty is the negative
#' log-density, at the given hyperparameters, with `NA` in the coordinates
#' that the distribution does not determine.
#'
#' @details
#' A proper penalty is a negative log-prior, and [penalty_value()] includes
#' its normalizing constant for that reason. Simulating from a model that
#' contains a random effect, a ridge or a lasso therefore amounts to drawing
#' those coefficients from this distribution.
#'
#' # Branches that draw
#'
#' A separable penalty draws from its own parent family, so a Gaussian prior
#' gives Gaussian effects, a Laplace prior gives Laplace effects and a Student
#' t prior gives heavy-tailed effects. A univariate parent is drawn once per
#' coordinate, and a multivariate parent once per block of coordinates.
#'
#' A quadratic, additive or structured penalty is a Gaussian prior with
#' precision \eqn{S(\theta)}, and two cases are drawn directly. Where \eqn{S}
#' is diagonal each coordinate is drawn on its own with standard deviation
#' \eqn{S_{jj}^{-1/2}}, which covers a ridge and a smooth in a basis whose
#' penalty is diagonal (such as the Demmler-Reinsch basis); a zero on that
#' diagonal is a direction that the prior does not determine, and its entry
#' is `NA`. Where \eqn{S} is not diagonal but has full rank, one Cholesky
#' factor \eqn{S = R^\top R} gives the draw, since \eqn{\beta = R^{-1}z} has
#' covariance \eqn{S^{-1}}.
#'
#' A rank-deficient \eqn{S} that is not diagonal, such as the penalty of an
#' anisotropic tensor smooth, gives a vector of `NA`. The prior is flat along
#' the null space of \eqn{S}, and a draw on the range alone would require a
#' basis of that range, a decomposition that these branches do not compute.
#'
#' SCAD and MCP are improper and are not densities, so their draws are `NA`
#' everywhere. A branch that is not quadratic and does not register a method
#' of its own also gives `NA` everywhere.
#'
#' # The map
#'
#' A separable penalty is a prior on \eqn{D\beta}, so a draw of the
#' coefficients inverts the map. Without a map the draw is returned as it is.
#' Under a diagonal map, which is what standardization produces, each
#' coordinate is divided by its diagonal entry, and a coordinate whose entry
#' is zero is not determined by the prior and is `NA`. Under any other map,
#' square and invertible included, the draw is not inverted and every
#' coordinate is `NA`. A quadratic penalty does not need this step, because
#' [penalty_matrix()] returns the precision in the coefficients themselves.
#'
#' # Undetermined coordinates
#'
#' The coordinates that the prior does not determine are `NA` rather than
#' zero. A caller that draws a whole model has its own rule for a coefficient
#' that is not covered by any prior, an intercept being the ordinary case.
#' A zero would be indistinguishable from a prior that concentrates at zero,
#' and the caller could not tell which coordinates remain to be filled.
#'
#' @param pen A [penalty()] object.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, as for the other generics of the package.
#' @param ... Passed to methods. The shipped methods do not read it.
#'
#' @return A numeric vector of length `pen@n_coef`, `NA` in the coordinates
#'   that the prior does not determine and possibly `NA` throughout.
#'
#' @examples
#' set.seed(1)
#'
#' # A Gaussian prior on ten effects gives ten Gaussian draws.
#' ridge <- ridge_penalty(n_coef = 10)
#' round(penalty_draw(ridge, list(lambda = 4)), 3)
#'
#' # The scale is that of the prior, so a larger lambda shrinks the draw.
#' c(loose = stats::sd(penalty_draw(ridge, list(lambda = 0.25))),
#'   tight = stats::sd(penalty_draw(ridge, list(lambda = 100))))
#'
#' # The lasso is a Laplace prior and gives Laplace draws.
#' round(penalty_draw(lasso_penalty(n_coef = 6), list(lambda = 2)), 3)
#'
#' # The unpenalized direction of a smooth is NA.
#' P <- diag(c(0, rep(1, 4)))
#' penalty_draw(quadratic_penalty(P), list(lambda = 1))
#'
#' # A diagonal map with a zero entry leaves that coordinate NA.
#' penalty_draw(lasso_penalty(map = Matrix::Diagonal(x = c(1, 0, 2))),
#'              list(lambda = 1))
#'
#' # SCAD is not a density, so every coordinate is NA.
#' penalty_draw(scad_penalty(n_coef = 3), list(lambda = 1, a = 3.7))
#'
#' @section Methods:
#' The method on the base class [penalty()] draws from a quadratic penalty
#' (for which [is_quadratic()] is `TRUE`) read as a Gaussian prior, as
#' described above, and returns `NA` for any other branch. A separable
#' penalty registers [penalty_draw.DistribPenalty()].
#'
#' @seealso [penalty_value()] for the negative log-density from which the
#'   draw is made, [is_proper()] for whether the penalty is a proper prior,
#'   [penalty_matrix()] for the precision of the Gaussian branches.
#'
#' @aliases penalty_draw.penalty
#' @export
penalty_draw <- S7::new_generic("penalty_draw", "pen",
  function(pen, theta, ...) {
    theta <- align_ptheta(pen, theta)
    S7::S7_dispatch()
  })

S7::method(penalty_draw, penalty) <- function(pen, theta, ...) {
  n <- as.integer(pen@n_coef)
  out <- rep(NA_real_, n)
  if (!is_quadratic(pen)) return(out)
  S <- penalty_matrix(pen, theta)
  # A diagonal precision is coordinatewise and is the ordinary case: a ridge,
  # and a Demmler-Reinsch smooth, whose penalty is diag(0, 1, ..., 1). The
  # zeros are the unpenalized directions and stay NA.
  if (isTRUE(suppressWarnings(try(Matrix::isDiagonal(S), silent = TRUE)))) {
    d <- as.numeric(Matrix::diag(S))
    pos <- is.finite(d) & d > 0
    if (any(pos)) out[pos] <- stats::rnorm(sum(pos), 0, 1 / sqrt(d[pos]))
    return(out)
  }
  # Otherwise only a proper prior can be drawn from without a basis of the
  # range, and there one triangular factor does it: with S = R'R the vector
  # R^-1 z has covariance S^-1.
  if (!is_proper(pen)) return(out)
  R <- try(chol(as.matrix(S)), silent = TRUE)
  if (inherits(R, "try-error")) return(out)
  as.numeric(backsolve(R, stats::rnorm(n)))
}

#' @title Drawing From a Separable Penalty
#' @name penalty_draw.DistribPenalty
#'
#' @description
#' Draws from the parent family, one draw per coordinate for a univariate
#' parent and one per block for a multivariate parent, carried back through
#' the map.
#'
#' @details
#' The penalty is the parent's negative log-density applied to each
#' coordinate (or block) of \eqn{D\beta}, so a draw is
#' [distributions7::distrib_rng()] at the hyperparameters, which are the
#' parent's own parameters. Under a multivariate parent the draw has one row
#' per block and is flattened in the order in which [dp_arg()] reads it.
#'
#' A diagonal map is inverted by dividing, since the prior is on
#' \eqn{d_j\beta_j}. A coordinate whose \eqn{d_j} is zero is not determined by
#' the prior and is `NA`, and the other coordinates are drawn as usual. A
#' square map whose off-diagonal entries are all zero is stored as a diagonal
#' map by [distrib_penalty()], whatever its class, and is treated in the same
#' way. Under any other map every coordinate is `NA`.
#'
#' @param pen A [DistribPenalty()] object.
#' @param theta Its hyperparameters, which are the parent's parameters.
#' @param ... Ignored.
#'
#' @return A numeric vector of length `pen@n_coef`, `NA` where a diagonal map
#'   has a zero entry and `NA` throughout under a map that is not diagonal.
#'
#' @seealso [penalty_draw()], [distrib_penalty()]
#' @keywords internal
S7::method(penalty_draw, DistribPenalty) <- function(pen, theta, ...) {
  n <- as.integer(pen@n_coef)
  d <- NULL
  if (!is.null(pen@map)) {
    d <- dp_diag_entries(pen@map)
    if (is.null(d)) return(rep(NA_real_, n))
  }
  t <- if (pen@block == 1L) {
    as.numeric(distributions7::distrib_rng(pen@parent, n, theta))
  } else {
    dp_flat(pen, distributions7::distrib_rng(pen@parent,
                                             n %/% pen@block, theta))
  }
  if (is.null(d)) return(t)
  # a zero on the diagonal is a coordinate the prior says nothing about
  out <- t / d
  out[is.na(d) | d == 0] <- NA_real_
  out
}
