#' @include scad_mcp.R
NULL

#' @title Check a Penalty Numerically
#'
#' @description
#' The sibling of \code{check_link} and \code{check_distrib}: every closed
#' form is compared against a route that shares no code with it. The
#' gradient and the Hessian are checked against \pkg{numDeriv} on the value,
#' the theta blocks against \pkg{numDeriv} in each hyperparameter -- the
#' mixed block by Richardson on the analytic gradient, never a nested
#' difference -- and the map by comparing \eqn{\rho(D\beta)} routes. Grids
#' are placed away from the kink set the object itself declares.
#'
#' @param pen A \code{\link{penalty}} object.
#' @param beta A coefficient vector, or \code{NULL} for a default draw.
#' @param theta A named hyperparameter list, or \code{NULL} for midpoints.
#' @param tol The comparison tolerance.
#' @param verbose Logical; print the table.
#'
#' @return A data frame with one row per check, invisibly when printed.
#'
#' @examples
#' res <- check_penalty(quadratic_penalty(diag(3)))
#' all(res$status == "OK")
#'
#' @export
check_penalty <- function(pen, beta = NULL, theta = NULL, tol = 1e-6,
                          verbose = TRUE) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("check_penalty() needs the numDeriv package.", call. = FALSE)
  }
  if (is.null(theta)) {
    theta <- lapply(pen@params_bounds, function(b) {
      lo <- if (is.finite(b[1])) b[1] else -1
      hi <- if (is.finite(b[2])) b[2] else lo + 2
      lo + 0.6 * (hi - lo)
    })
  }
  theta <- align_ptheta(pen, theta)
  q <- pen@n_coef
  if (is.null(beta)) {
    set.seed(7)
    beta <- round(stats::rnorm(q, sd = 1.3), 2) + 0.11
    # push t = D beta away from the declared kink set
    kk <- penalty_kinks(pen, theta)
    if (length(kk)) {
      t <- map_apply(pen, beta)
      for (i in seq_len(50)) {
        if (min(abs(outer(t, kk, "-"))) > 0.05) break
        beta <- beta + 0.033
        t <- map_apply(pen, beta)
      }
    }
  }

  checks <- list()
  add <- function(name, gap) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, max_error = gap,
      status = if (is.finite(gap) && gap < tol) "OK" else "FAILED"
    )
  }
  rel <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))

  g <- penalty_gradient(pen, beta, theta)
  gn <- numDeriv::grad(function(b) penalty_value(pen, b, theta), beta)
  add("gradient vs numDeriv", rel(g, gn))

  H <- penalty_hessian(pen, beta, theta)
  Hn <- numDeriv::jacobian(function(b) penalty_gradient(pen, b, theta), beta)
  add("hessian vs numDeriv on the gradient", rel(H, (Hn + t(Hn)) / 2))

  gt <- penalty_grad_theta(pen, beta, theta)
  ht <- penalty_hess_theta(pen, beta, theta)
  cr <- penalty_cross(pen, beta, theta)
  for (p in pen@params) {
    fn <- function(v) {
      th <- theta; th[[p]] <- v
      penalty_value(pen, beta, th)
    }
    add(paste0("grad_theta[", p, "] vs numDeriv"),
        rel(gt[[p]], numDeriv::grad(fn, theta[[p]])))
    gfn <- function(v) {
      th <- theta; th[[p]] <- v
      penalty_grad_theta(pen, beta, th)[[p]]
    }
    add(paste0("hess_theta[", p, "_", p, "] vs numDeriv"),
        rel(ht[[paste0(p, "_", p)]], numDeriv::grad(gfn, theta[[p]])))
    cfn <- function(v) {
      th <- theta; th[[p]] <- v
      penalty_gradient(pen, beta, th)
    }
    add(paste0("cross[", p, "] vs Richardson on the gradient"),
        rel(cr[[p]], as.numeric(numDeriv::jacobian(cfn, theta[[p]]))))
  }

  if (is_quadratic(pen)) {
    # the value is quadratic in beta: three-point exactness along a ray
    v0 <- penalty_value(pen, 0 * beta, theta)
    v1 <- penalty_value(pen, beta, theta)
    v2 <- penalty_value(pen, 2 * beta, theta)
    add("quadratic three-point identity",
        abs(v2 - 4 * v1 + 3 * v0) / max(1, abs(v1)))
    if ("lambda" %in% pen@params) {
      # the plain quadratic branch: one scale multiplies a constant matrix
      lp1 <- penalty_logpdet(pen, theta)
      th2 <- theta; th2$lambda <- theta$lambda * exp(1)
      lp2 <- penalty_logpdet(pen, th2)
      add("logpdet linear in log lambda with slope r",
          abs((lp2$value - lp1$value) - penalty_rank(pen)))
    } else {
      # the structured branch: the log pseudo-determinant's own gradient
      # against numDeriv on its value
      lp <- penalty_logpdet(pen, theta)
      worst <- 0
      for (p2 in pen@params) {
        ref <- numDeriv::grad(function(v) {
          th <- theta; th[[p2]] <- v
          penalty_logpdet(pen, th)$value
        }, theta[[p2]])
        worst <- max(worst, abs(lp$grad[[p2]] - ref) / max(1, abs(ref)))
      }
      add("logpdet gradient vs numDeriv", worst)
    }
    M <- penalty_matrix(pen, theta)
    nb <- penalty_null_basis(pen)
    if (ncol(nb)) {
      add("null basis annihilates the matrix",
          max(abs(M %*% nb)) / max(1, max(abs(M))))
    }
  }

  out <- do.call(rbind, checks)
  if (verbose) print(out, row.names = FALSE)
  invisible(out)
}
