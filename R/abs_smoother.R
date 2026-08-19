#' @include penalty_class.R
NULL

#' @title S7 Class for Smoothers of the Absolute Value
#'
#' @description
#' A smooth replacement \eqn{s(u)} for \eqn{\lvert u \rvert}, carrying its
#' derivatives in \eqn{u} up to order five as functions. Every non-smooth
#' primitive the toolkit uses is generated from the absolute value --
#'
#' \deqn{\operatorname{sign}(u) = \frac{d\lvert u\rvert}{du}, \qquad
#'   \mathbb{1}(u \geq 0) = \frac{1 + \operatorname{sign}(u)}{2}, \qquad
#'   (u)_{+} = \frac{u + \lvert u\rvert}{2},}
#'
#' so one contract serves them all: the smooth sign is \eqn{s'}, the smooth
#' step \eqn{(1 + s'(u))/2} and the smooth hinge \eqn{(u + s(u))/2} follow by
#' composition, and a consumer that replaces \eqn{\lvert u\rvert} by
#' \eqn{s(u)} has replaced every one of them consistently. A break-point term
#' smoothed this way becomes an ordinary nonlinear term whose block is the
#' true Jacobian, which is what makes a random or penalized development of
#' its break-points fittable; a kinked penalty smoothed this way becomes a
#' proper separable one, at the price of the exact zeros the kink produced.
#'
#' @details
#' The derivatives are FUNCTIONS and not expressions: a piecewise smoother
#' (the quintic) has branches, which \code{\link[stats]{deriv}} does not
#' read, and with the derivatives written in the constructor the branches
#' are ordinary code. Each function takes \code{(u, width)} and vectorizes
#' in both, so a per-group width is one value per observation.
#'
#' \code{width} is the transition scale -- \code{h} for a smoother whose
#' parameter is a length, the bent-cable reading of a transition of width
#' \eqn{h}; \code{c} for the hyperbolic, whose parameter is a squared
#' length. \code{NULL}, the default, asks the consumer to resolve it from
#' the data at build (a break-point term takes the median spacing of its
#' covariate), through \code{\link{smoother_width}}. \code{per_group} asks
#' for one width per group where a grouping is available, the validity
#' window of a Laplace approximation being per-subject.
#'
#' \code{tau_correction} is a property of the mollifier, not of any model:
#' the probit smoother satisfies an exact convolution identity, smoothing
#' with width \eqn{h} being the same as convolving the break-point with
#' \eqn{N(0, h^2)}, so an apparent scale \eqn{\tau} of a random break-point
#' composes as \eqn{\tau^2_{\mathrm{true}} = \tau^2 - h^2} and the smoother
#' declares the correction. Smoothers with no such identity declare
#' \code{NULL} and a consumer reports the apparent scale alone.
#'
#' @param smoother_name A string naming the smoother.
#' @param width The transition width, or \code{NULL} to be resolved at
#'   build.
#' @param width_name What the width parameter is called (\code{"h"} or
#'   \code{"c"}).
#' @param per_group Whether the width is resolved per group.
#' @param s A list of six functions of \code{(u, width)}: \eqn{s} and its
#'   derivatives in \eqn{u} of orders one to five, in order.
#' @param width_from_spacing \code{NULL} (the identity) or a function
#'   carrying a spacing in covariate units onto the width parameter's own
#'   scale (the square, for the hyperbolic).
#' @param tau_correction \code{NULL}, or a function \code{(tau, width)}
#'   returning the corrected scale of a random break-point.
#' @param exact_radius \code{NULL}, or a function of the width returning
#'   the radius beyond which \eqn{s(u) = \lvert u\rvert} exactly.
#'
#' @return An object of class \code{abs_smoother}.
#'
#' @references
#' Bacon, D. W. and Watts, D. G. (1971). Estimating the transition between
#' two intersecting straight lines. \emph{Biometrika}, 58(3), 525--534.
#'
#' Tishler, A. and Zang, I. (1981). A new maximum likelihood algorithm for
#' piecewise regression. \emph{Journal of the American Statistical
#' Association}, 76(376), 980--987.
#'
#' Seo, M. H. and Linton, O. (2007). A smoothed least squares estimator
#' for threshold regression models. \emph{Journal of Econometrics},
#' 141(2), 704--735.
#'
#' @seealso \code{\link{smooth_probit}}, \code{\link{smooth_hyperbolic}},
#'   \code{\link{smooth_quintic}}, \code{\link{check_abs_smoother}}
#'
#' @examples
#' S7::S7_inherits(smooth_probit(), abs_smoother)
#'
#' @export
abs_smoother <- S7::new_class(
  name = "abs_smoother",
  properties = list(
    smoother_name = S7::class_character,
    width = S7::class_any,
    width_name = S7::new_property(S7::class_character, default = "h"),
    per_group = S7::new_property(S7::class_logical, default = FALSE),
    s = S7::class_list,
    width_from_spacing = S7::class_any,
    tau_correction = S7::class_any,
    exact_radius = S7::class_any
  ),
  validator = function(self) {
    if (length(self@smoother_name) != 1L || is.na(self@smoother_name) ||
        !nzchar(self@smoother_name)) {
      return("@smoother_name must be a single non-empty string")
    }
    if (!is.null(self@width) &&
        (!is.numeric(self@width) || length(self@width) != 1L ||
         !is.finite(self@width) || self@width <= 0)) {
      return("@width must be NULL or a single positive number")
    }
    if (length(self@width_name) != 1L || !nzchar(self@width_name)) {
      return("@width_name must be a single non-empty string")
    }
    if (length(self@per_group) != 1L || is.na(self@per_group)) {
      return("@per_group must be TRUE or FALSE")
    }
    if (length(self@s) != 6L ||
        !all(vapply(self@s, is.function, logical(1)))) {
      return(paste("@s must be a list of six functions of (u, width): s and",
                   "its derivatives of orders one to five"))
    }
    if (any(vapply(self@s, function(f) length(formals(f)) < 2L,
                   logical(1)))) {
      return("every function in @s must take (u, width)")
    }
    for (nm in c("width_from_spacing", "tau_correction", "exact_radius")) {
      v <- S7::prop(self, nm)
      if (!is.null(v) && !is.function(v)) {
        return(sprintf("@%s must be NULL or a function", nm))
      }
    }
    NULL
  }
)

#' The Probit Smoother of the Absolute Value
#'
#' @description
#' \eqn{s(u) = \mathbb{E}\lvert u + hZ\rvert} for \eqn{Z} standard normal:
#'
#' \deqn{s(u) = u\,\bigl(2\Phi(u/h) - 1\bigr) + 2h\,\phi(u/h), \qquad
#'   s'(u) = 2\Phi(u/h) - 1, \qquad s''(u) = \frac{2}{h}\,\phi(u/h),}
#'
#' the recommended default. Its excess over \eqn{\lvert u\rvert} has
#' gaussian tails, so the smoothing bias is confined to a window of width
#' \eqn{h} around the kink, and it satisfies an exact convolution identity:
#' smoothing the step with width \eqn{h} is convolving the break-point with
#' \eqn{N(0, h^2)}, so the apparent scale of a random break-point composes
#' as \eqn{\tau^2_{\mathrm{apparent}} = \tau^2 + h^2} and the smoother
#' declares the closed correction
#' \eqn{\tau_{\mathrm{true}} = \sqrt{\tau^2 - h^2}}.
#'
#' @param h The transition width, or \code{NULL} (the default) to be
#'   resolved at build from the covariate's spacing.
#' @param per_group Whether the width is resolved per group.
#'
#' @return An \code{\link{abs_smoother}}.
#'
#' @examples
#' sm <- smooth_probit(h = 0.2)
#' smoother_deriv(sm, 0, order = 0)  # 2 * h * dnorm(0)
#'
#' @seealso \code{\link{abs_smoother}}, \code{\link{smooth_hyperbolic}},
#'   \code{\link{smooth_quintic}}
#' @export
smooth_probit <- function(h = NULL, per_group = FALSE) {
  abs_smoother(
    smoother_name = "probit", width = h, width_name = "h",
    per_group = per_group,
    s = list(
      function(u, width) {
        z <- u / width
        u * (2 * stats::pnorm(z) - 1) + 2 * width * stats::dnorm(z)
      },
      function(u, width) 2 * stats::pnorm(u / width) - 1,
      function(u, width) 2 * stats::dnorm(u / width) / width,
      function(u, width) {
        z <- u / width
        -2 * z * stats::dnorm(z) / width^2
      },
      function(u, width) {
        z <- u / width
        2 * (z^2 - 1) * stats::dnorm(z) / width^3
      },
      function(u, width) {
        z <- u / width
        2 * (3 * z - z^3) * stats::dnorm(z) / width^4
      }
    ),
    width_from_spacing = NULL,
    tau_correction = function(tau, width) sqrt(pmax(tau^2 - width^2, 0)),
    exact_radius = NULL
  )
}

#' The Hyperbolic Smoother of the Absolute Value
#'
#' @description
#' \eqn{s(u) = \sqrt{u^2 + c}}, the simplest smoother, whose parameter
#' \eqn{c} is a squared length (the transition width in covariate units is
#' \eqn{\sqrt{c}}). Its excess over \eqn{\lvert u\rvert} decays only as
#' \eqn{c/(4\lvert u\rvert)} in the tails, so the smoothing bias spreads
#' away from the kink -- the worst bias profile of the three at equal
#' width -- and no convolution identity relates \eqn{c} to the scale of a
#' random break-point, so \code{tau_correction} is \code{NULL} and a
#' consumer reports the apparent scale alone.
#'
#' @param c The squared transition width, or \code{NULL} (the default) to
#'   be resolved at build as the square of the covariate's spacing.
#'
#' @return An \code{\link{abs_smoother}}.
#'
#' @examples
#' sm <- smooth_hyperbolic(c = 0.04)
#' smoother_deriv(sm, 0, order = 0)  # sqrt(c)
#'
#' @seealso \code{\link{abs_smoother}}, \code{\link{smooth_probit}}
#' @export
smooth_hyperbolic <- function(c = NULL) {
  abs_smoother(
    smoother_name = "hyperbolic", width = c, width_name = "c",
    per_group = FALSE,
    s = list(
      function(u, width) sqrt(u^2 + width),
      function(u, width) u / sqrt(u^2 + width),
      function(u, width) width / (u^2 + width)^1.5,
      function(u, width) -3 * width * u / (u^2 + width)^2.5,
      function(u, width) 3 * width * (4 * u^2 - width) / (u^2 + width)^3.5,
      function(u, width) {
        15 * width * u * (3 * width - 4 * u^2) / (u^2 + width)^4.5
      }
    ),
    width_from_spacing = function(gap) gap^2,
    tau_correction = NULL,
    exact_radius = NULL
  )
}

#' The Quintic Smoother of the Absolute Value
#'
#' @description
#' \eqn{\lvert u\rvert} replaced inside \eqn{[-h, h]} by the even
#' polynomial matching it to third order at the seam, and EXACT outside:
#' with \eqn{t = u/h},
#'
#' \deqn{s(u) = \frac{h}{16}\,\bigl(5 + 15t^2 - 5t^4 + t^6\bigr)
#'   \quad \text{for } \lvert t\rvert < 1, \qquad
#'   s(u) = \lvert u\rvert \quad \text{otherwise},}
#'
#' so \eqn{s''(u) = 15(1 - t^2)^2/(8h)} vanishes to second order at
#' \eqn{\pm h} and the function is \eqn{C^3}, with a jump in the fourth
#' derivative at the seam. The smoothing bias is exactly zero outside
#' \eqn{[\psi - h, \psi + h]}, which is the cleanest fixed-width choice;
#' the branches are why the contract carries the derivatives as functions,
#' \code{\link[stats]{deriv}} not reading a clamp.
#'
#' @param h The transition half-width, or \code{NULL} (the default) to be
#'   resolved at build from the covariate's spacing.
#'
#' @return An \code{\link{abs_smoother}}.
#'
#' @examples
#' sm <- smooth_quintic(h = 0.5)
#' smoother_deriv(sm, 1, order = 0)  # exactly 1: outside the transition
#'
#' @seealso \code{\link{abs_smoother}}, \code{\link{smooth_probit}}
#' @export
smooth_quintic <- function(h = NULL) {
  inside <- function(u, width) abs(u) < width
  abs_smoother(
    smoother_name = "quintic", width = h, width_name = "h",
    per_group = FALSE,
    s = list(
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width),
               width * (5 + 15 * t^2 - 5 * t^4 + t^6) / 16, abs(u))
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), (15 * t - 10 * t^3 + 3 * t^5) / 8, sign(u))
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), 15 * (1 - t^2)^2 / (8 * width), 0)
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), -15 * t * (1 - t^2) / (2 * width^2), 0)
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), 15 * (3 * t^2 - 1) / (2 * width^3), 0)
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), 45 * t / width^4, 0)
      }
    ),
    width_from_spacing = NULL,
    tau_correction = NULL,
    exact_radius = function(width) width
  )
}

#' Evaluate a Smoother or One of Its Derivatives
#'
#' @description
#' \eqn{d^k s/du^k} at the given points, for \eqn{k} from 0 to 5. The width
#' is the smoother's own unless one is supplied, which is how a consumer
#' that resolved the width at build evaluates without mutating the object;
#' a vector width is one value per point.
#'
#' @param smoother An \code{\link{abs_smoother}}.
#' @param u A numeric vector.
#' @param width The width, or \code{NULL} for the smoother's own.
#' @param order 0 to 5.
#'
#' @return A numeric vector as long as \code{u}.
#'
#' @examples
#' smoother_deriv(smooth_hyperbolic(c = 1), 0:3, order = 1)
#'
#' @seealso \code{\link{abs_smoother}}
#' @export
smoother_deriv <- function(smoother, u, width = NULL, order = 0L) {
  if (!S7::S7_inherits(smoother, abs_smoother)) {
    stop("'smoother' must be an abs_smoother.", call. = FALSE)
  }
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 0L || order > 5L) {
    stop("'order' must be an integer from 0 to 5.", call. = FALSE)
  }
  if (is.null(width)) width <- smoother@width
  if (is.null(width)) {
    stop(paste("the smoother's width is unresolved; supply 'width' or",
               "construct the smoother with one."), call. = FALSE)
  }
  smoother@s[[order + 1L]](as.numeric(u), width)
}

#' Resolve a Smoother's Width from a Spacing
#'
#' @description
#' The width the smoother carries where it carries one, and otherwise the
#' given spacing carried onto the width parameter's own scale (the square,
#' for a smoother parametrized by a squared length). It is what a consumer
#' calls at build: a break-point term hands it the median spacing of its
#' covariate, the smallest transition the data can tell from a step.
#'
#' @param smoother An \code{\link{abs_smoother}}.
#' @param spacing A spacing in covariate units, one value or one per group.
#'
#' @return The width, the same length as \code{spacing} when resolved from
#'   it.
#'
#' @examples
#' smoother_width(smooth_probit(), 0.3)          # 0.3
#' smoother_width(smooth_hyperbolic(), 0.3)      # 0.09
#' smoother_width(smooth_probit(h = 0.5), 0.3)   # 0.5, the width it holds
#'
#' @seealso \code{\link{smoother_width_floor}}
#' @export
smoother_width <- function(smoother, spacing) {
  if (!S7::S7_inherits(smoother, abs_smoother)) {
    stop("'smoother' must be an abs_smoother.", call. = FALSE)
  }
  if (!is.null(smoother@width)) return(smoother@width)
  spacing <- as.numeric(spacing)
  if (!length(spacing) || anyNA(spacing) || any(spacing <= 0)) {
    stop("'spacing' must be positive.", call. = FALSE)
  }
  f <- smoother@width_from_spacing
  if (is.null(f)) spacing else f(spacing)
}

#' The Smallest Width a Consumer May Use
#'
#' @description
#' The floor of the width, derived from the expression that binds rather
#' than chosen. The derivatives of a smoother scale as
#' \eqn{s^{(k)} \sim h^{1-k}}, so the Jacobian column of a smoothed
#' break-point carries \eqn{s''(0)/2 \sim 1/h} against covariate columns of
#' order the range \eqn{D}: holding the design's condition number below
#' \eqn{\epsilon^{-1/2}}, which leaves half the digits of a double to a QR
#' of the design -- the same bound the working schedule's scaling floor
#' rests on -- gives \eqn{h \geq \sqrt{\epsilon}\,D}, carried onto the
#' width parameter's own scale for a smoother parametrized by a squared
#' length.
#'
#' @param smoother An \code{\link{abs_smoother}}.
#' @param scale The scale of the covariate (its range).
#'
#' @return A single number.
#'
#' @examples
#' smoother_width_floor(smooth_probit(), scale = 10)
#'
#' @seealso \code{\link{smoother_width}}
#' @export
smoother_width_floor <- function(smoother, scale) {
  if (!S7::S7_inherits(smoother, abs_smoother)) {
    stop("'smoother' must be an abs_smoother.", call. = FALSE)
  }
  scale <- as.numeric(scale)
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("'scale' must be a single positive number.", call. = FALSE)
  }
  h <- sqrt(.Machine$double.eps) * scale
  f <- smoother@width_from_spacing
  if (is.null(f)) h else f(h)
}

S7::method(print, abs_smoother) <- function(x, ...) {
  cat(sprintf("<abs_smoother> %s: %s%s\n", x@smoother_name,
              if (is.null(x@width)) {
                sprintf("%s resolved at build", x@width_name)
              } else {
                sprintf("%s = %g", x@width_name, x@width)
              },
              if (isTRUE(x@per_group)) ", per group" else ""))
  if (!is.null(x@tau_correction)) {
    cat("  declares a scale correction for a random break-point\n")
  }
  if (!is.null(x@exact_radius)) {
    cat("  exact outside the transition\n")
  }
  invisible(x)
}

#' Check a Smoother of the Absolute Value Numerically
#'
#' @description
#' The sibling of \code{check_penalty} for user-written smoothers: the
#' structural properties of \eqn{s} (even, with an odd first derivative
#' bounded by one, convex, and matching \eqn{\lvert u\rvert} in the tails)
#' and each derivative order against one numerical differentiation of the
#' analytic order below it, never a nested difference.
#'
#' @details
#' The grid covers the transition and the tails, placed by the smoother's
#' own intrinsic scale \eqn{s(0)}. Order \eqn{k} is compared against
#' \pkg{numDeriv} applied ONCE to the analytic order \eqn{k - 1}, so a
#' wrong derivative is caught against an independent route while the
#' reference never degenerates into a difference of differences.
#'
#' @param smoother An \code{\link{abs_smoother}}.
#' @param width The width to check at, or \code{NULL} for the smoother's
#'   own (0.5 where it carries none).
#' @param tol The comparison tolerance.
#' @param verbose Logical; print the table.
#'
#' @return A data frame with one row per check, invisibly when printed.
#'
#' @examples
#' res <- check_abs_smoother(smooth_probit(h = 0.3))
#' all(res$status == "OK")
#'
#' @seealso \code{\link{abs_smoother}}, \code{\link{check_penalty}}
#' @export
check_abs_smoother <- function(smoother, width = NULL, tol = 1e-6,
                               verbose = TRUE) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("check_abs_smoother() needs the numDeriv package.", call. = FALSE)
  }
  if (!S7::S7_inherits(smoother, abs_smoother)) {
    stop("'smoother' must be an abs_smoother.", call. = FALSE)
  }
  if (is.null(width)) width <- smoother@width
  if (is.null(width)) width <- 0.5

  # the intrinsic scale: s(0) is of order the transition width for any
  # smoother of the absolute value, so the grid needs no reading of what
  # the width parameter means
  w0 <- max(smoother_deriv(smoother, 0, width, 0L), sqrt(width) * 1e-3)
  # off the seam and off zero: a piecewise smoother has measure-zero points
  # where a one-sided derivative is read, and a difference straddling one
  # compares nothing
  grid <- w0 * c(-6.3, -2.7, -1.13, -0.41, 0.17, 0.59, 1.21, 3.1, 7.7)

  checks <- list()
  add <- function(name, gap) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, max_error = gap,
      status = if (is.finite(gap) && gap < tol) "OK" else "FAILED"
    )
  }

  s0 <- smoother_deriv(smoother, grid, width, 0L)
  add("s is even", max(abs(s0 - smoother_deriv(smoother, -grid, width, 0L))) /
        max(1, max(abs(s0))))
  s1 <- smoother_deriv(smoother, grid, width, 1L)
  add("s' is odd", max(abs(s1 + smoother_deriv(smoother, -grid, width, 1L))))
  add("|s'| bounded by one", max(0, max(abs(s1)) - 1))
  add("s convex", max(0, -min(smoother_deriv(smoother, grid, width, 2L))))
  # the excess over |u| leaves the transition: at ten intrinsic widths it is
  # at most a fifth of what it is at the kink, which the hyperbolic's
  # polynomial tails still satisfy (excess c / (4|u|), a twentieth there)
  e0 <- smoother_deriv(smoother, 0, width, 0L)
  e10 <- smoother_deriv(smoother, 10 * w0, width, 0L) - 10 * w0
  add("matches |u| in the tails", max(0, e10 - e0 / 5))

  for (k in 1:5) {
    fk <- vapply(grid, function(u) {
      numDeriv::grad(function(v) smoother_deriv(smoother, v, width, k - 1L),
                     u)
    }, numeric(1))
    ak <- smoother_deriv(smoother, grid, width, k)
    add(sprintf("order %d vs numDeriv on order %d", k, k - 1L),
        max(abs(ak - fk)) / max(1, max(abs(fk))))
  }

  tc <- smoother@tau_correction
  if (!is.null(tc)) {
    tau <- 3 * w0
    add("tau_correction bounded by tau", max(0, tc(tau, width) - tau))
    add("tau_correction is the identity at width zero",
        abs(tc(tau, width * 1e-8) - tau) / tau)
  }

  out <- do.call(rbind, checks)
  rownames(out) <- NULL
  if (verbose) {
    cat(sprintf("check_abs_smoother: %s at %s = %g\n",
                smoother@smoother_name, smoother@width_name, width))
    print(out, row.names = FALSE)
    return(invisible(out))
  }
  out
}
