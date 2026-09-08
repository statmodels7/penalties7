#' @include generics.R distrib_penalty.R quadratic_penalty.R prox.R
NULL

#' @title Drawing From a Penalty Read as a Prior
#'
#' @description
#' A vector drawn from the distribution whose negative log-density the penalty
#' is, at the hyperparameters given, with `NA` in the coordinates that
#' distribution does not determine.
#'
#' @details
#' A penalty is a negative log-prior, which is why [penalty_value()] keeps the
#' normalizing constant, and a prior is something one draws from. Simulating
#' from a model that carries a random effect, a ridge or a lasso means drawing
#' those coefficients from exactly this distribution instead of from a normal
#' chosen by whoever wrote the simulation.
#'
#' # Which branches answer
#'
#' A separable penalty draws coordinatewise from its own parent family, so a
#' Gaussian prior gives Gaussian effects, a Laplace prior gives Laplace ones
#' and a Student t prior gives heavy-tailed ones. That is the case a random
#' effect is, and the one this exists for.
#'
#' A quadratic, additive or structured penalty is a Gaussian prior with
#' precision \eqn{S(\theta)}, and two of its shapes are cheap. Where \eqn{S} is
#' diagonal each coordinate is drawn on its own at a standard deviation
#' \eqn{S_{jj}^{-1/2}}, which covers a ridge and a Demmler-Reinsch smooth,
#' whose penalty is \eqn{\mathrm{diag}(0, 1, \dots, 1)}; a zero on that
#' diagonal is a direction the prior says nothing about and comes back `NA`.
#' Where \eqn{S} is not diagonal but has full rank, one Cholesky factor gives
#' the draw, \eqn{\beta = R^{-1}z} having covariance \eqn{S^{-1}}.
#'
#' A deficient non-diagonal \eqn{S}, which an anisotropic tensor smooth has,
#' comes back all `NA`. The prior is flat along its null space and drawing on
#' the range alone needs a basis of that range, which is a decomposition of a
#' matrix these branches are built not to form.
#'
#' SCAD and MCP are improper by construction, densities of nothing, and answer
#' `NA` everywhere. So does any penalty whose branch says nothing here.
#'
#' # The map
#'
#' The penalty is a prior on \eqn{D\beta}, so a draw of the coefficients
#' inverts the map. For a separable penalty the identity and a diagonal map,
#' which is what standardization is, invert coordinatewise, and anything else
#' comes back `NA`: under a general map the prior leaves the coefficients
#' underdetermined, which is the same fact that makes [penalty_prox()] reject
#' one. A quadratic penalty needs no such care, [penalty_matrix()] returning
#' the precision in the coefficients themselves.
#'
#' # Why the unanswered coordinates are `NA` rather than zero
#'
#' A caller drawing a whole model has a rule of its own for a coefficient no
#' prior covers, an intercept being the ordinary case. Returning zero would be
#' indistinguishable from a prior that genuinely concentrates there, and the
#' caller could not tell which coordinates it still has to fill.
#'
#' @param pen A [penalty()] object.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same, as every other generic here takes them.
#' @param ... Passed to methods. No shipped method reads it.
#'
#' @return A numeric vector of length `pen@n_coef`, `NA` in the coordinates
#'   the prior does not determine and possibly `NA` throughout.
#'
#' @examples
#' set.seed(1)
#'
#' # A Gaussian prior on ten effects: ten Gaussian effects.
#' ridge <- ridge_penalty(n_coef = 10)
#' round(penalty_draw(ridge, list(lambda = 4)), 3)
#'
#' # The scale is the prior's own, so a larger lambda shrinks the draw.
#' c(loose = stats::sd(penalty_draw(ridge, list(lambda = 0.25))),
#'   tight = stats::sd(penalty_draw(ridge, list(lambda = 100))))
#'
#' # A Laplace prior gives Laplace effects, which is what a lasso says.
#' round(penalty_draw(lasso_penalty(n_coef = 6), list(lambda = 2)), 3)
#'
#' # A smooth leaves its null space undetermined, and says so.
#' P <- diag(c(0, rep(1, 4)))
#' penalty_draw(quadratic_penalty(P), list(lambda = 1))
#'
#' # SCAD is the density of nothing.
#' penalty_draw(scad_penalty(n_coef = 3), list(lambda = 1, a = 3.7))
#'
#' @seealso [penalty_value()] for the log-density this inverts,
#'   [is_proper()] for whether there is a distribution to draw from at all,
#'   [penalty_matrix()] for the precision the Gaussian branches read.
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
#' One draw per coordinate from the parent family, carried back through the
#' map.
#'
#' @details
#' The penalty is the parent's negative log-density applied to each coordinate
#' of \eqn{D\beta}, so a draw is [distributions7::distrib_rng()] at the
#' hyperparameters, which are the parent's own parameters. Under a
#' multivariate parent the draw is one row per block, flattened the way
#' [dp_arg()] reads it.
#'
#' A diagonal map is inverted by dividing, the prior being on \eqn{d_j\beta_j};
#' any other map comes back `NA`, the coefficients not being determined by a
#' prior on a lower-dimensional image of them.
#'
#' @param pen A [DistribPenalty()] object.
#' @param theta Its hyperparameters, which are the parent's parameters.
#' @param ... Ignored.
#'
#' @return A numeric vector of length `pen@n_coef`, all `NA` under a
#'   non-diagonal map.
#'
#' @seealso [penalty_draw()], [distrib_penalty()]
#' @keywords internal
S7::method(penalty_draw, DistribPenalty) <- function(pen, theta, ...) {
  n <- as.integer(pen@n_coef)
  d <- NULL
  if (!is.null(pen@map)) {
    d <- map_diagonal(pen)
    if (is.null(d)) return(rep(NA_real_, n))
  }
  t <- if (pen@block == 1L) {
    as.numeric(distributions7::distrib_rng(pen@parent, n, theta))
  } else {
    dp_flat(pen, distributions7::distrib_rng(pen@parent,
                                             n %/% pen@block, theta))
  }
  if (is.null(d)) t else t / d
}
