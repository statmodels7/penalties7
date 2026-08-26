#' @include scad_mcp.R
NULL

#' @title Check a Penalty Numerically
#'
#' @description
#' Compares every closed form a penalty declares against a route that shares
#' no code with it, and returns one row per comparison carrying the worst
#' relative error and a pass or fail. The gradient and the Hessian go against
#' \pkg{numDeriv} on the value; each hyperparameter block goes against
#' \pkg{numDeriv} in that hyperparameter, the mixed block by Richardson on the
#' analytic gradient so that no difference is ever taken of another difference;
#' and a quadratic penalty has its three-point identity, its log
#' pseudo-determinant and its null basis tested as well. Write a penalty of
#' your own and this says whether its derivatives are right.
#'
#' @details
#' # The rows
#'
#' Two rows are produced for every penalty, three more for each
#' hyperparameter, and up to three more when [is_quadratic()] is `TRUE`:
#'
#' | row | compares |
#' |---|---|
#' | `gradient vs numDeriv` | [penalty_gradient()] against a numerical gradient of [penalty_value()] |
#' | `hessian vs numDeriv on the gradient` | [penalty_hessian()] against a numerical Jacobian of the analytic gradient, symmetrized |
#' | `grad_theta[p] vs numDeriv` | [penalty_grad_theta()] against a numerical derivative of the value in hyperparameter `p` |
#' | `hess_theta[p_p] vs numDeriv` | [penalty_hess_theta()] against a numerical derivative of the analytic `grad_theta` |
#' | `cross[p] vs Richardson on the gradient` | [penalty_cross()] against a numerical Jacobian of the analytic coefficient gradient in `p` |
#' | `quadratic three-point identity` | \eqn{\rho(2\beta) - 4\rho(\beta) + 3\rho(0)}, which vanishes when the value is a quadratic form with no linear term |
#' | `logpdet linear in log lambda with slope r` | raising `lambda` by a factor of \eqn{e} raises [penalty_logpdet()] by exactly the rank |
#' | `logpdet gradient vs numDeriv` | the `grad` element of `penalty_logpdet()` against a numerical derivative of its `value`, for a quadratic penalty whose hyperparameters are not a single `lambda` |
#' | `null basis annihilates the matrix` | \eqn{PN}, where `N` is [penalty_null_basis()] |
#'
#' The last two of the log-determinant rows are alternatives: the first is
#' taken when `lambda` is among the hyperparameters, which is the plain
#' quadratic branch of one scale multiplying a constant matrix, and the second
#' otherwise, which is the structured branch. The null-basis row appears only
#' when the null basis has columns.
#'
#' So the count is `2 + 3 * length(pen@params)` plus one, two or three. A
#' lasso gives 5 rows, a full-rank ridge 7, a quadratic penalty over second
#' differences 8, and a structured penalty over a `3 x 3` log-Cholesky
#' precision 22.
#'
#' # What is not checked
#'
#' The pass is over the value and its derivatives, and over the pieces a
#' marginal criterion reads. It does not touch [penalty_prox()],
#' [penalty_prox_spec()], [penalty_dhessian()], [penalty_d2hessian()],
#' [penalty_dcross()] or [penalty_readable()]. It calls [penalty_kinks()] to
#' place the grid but never tests the kinks themselves.
#'
#' An [additive_penalty()] reports `is_quadratic()` as `FALSE`, so none of the
#' three quadratic rows runs for it and its rank, matrix and log
#' pseudo-determinant go untested. A two-component additive penalty therefore
#' produces the bare 8 rows its two hyperparameters earn.
#'
#' # Where the grid is placed
#'
#' A penalty with a kink has no derivative there, so a numerical reference
#' straddling one measures the kink and not the formula. With `beta` left at
#' `NULL` the draw is nudged in steps of `0.033`, up to fifty times, until
#' every coordinate of \eqn{D\beta} sits at least `0.05` from every kink the
#' object declares. At the default coefficient count the first draw already
#' clears the SCAD and MCP kink sets and no step is taken.
#'
#' @param pen A [penalty()] object, of any branch.
#' @param beta A numeric coefficient vector of length `pen@n_coef`. `NULL`,
#'   the default, draws one from `rnorm(sd = 1.3)` rounded to two places and
#'   shifted by `0.11`, then pushes it clear of the kinks as above.
#'   The draw is taken from a fixed seed, so the report is the same on two runs,
#'   and the caller's own `.Random.seed` is restored on exit, so a call in the
#'   middle of a simulation leaves that simulation unchanged.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same. `NULL`, the default, places each hyperparameter
#'   six tenths of the way across its own bounds, reading an infinite lower
#'   bound as `-1` and an infinite upper bound as two above the lower. That
#'   gives `lambda = 1.2` on every branch that has one, `alpha = 0.6` for the
#'   elastic net, `a = 3.2` for SCAD and `gamma = 2.2` for MCP.
#' @param tol The relative error above which a row is reported as `FAILED`.
#'   A single positive number, `1e-6` by default. Over the nine shipped
#'   branches the worst error measured is `2.1e-10`, four orders under it, so
#'   the default separates a correct penalty from one whose formula is wrong
#'   in its fifth digit.
#' @param verbose `TRUE`, the default, prints the table without row names.
#'   The result is returned invisibly either way.
#'
#' @return A data frame with one row per check and three columns: `check`
#'   (character, the row's name as tabulated above), `max_error` (numeric, the
#'   worst absolute difference divided by `max(1, max(abs(reference)))`), and
#'   `status` (character, `"OK"` or `"FAILED"`). Returned invisibly.
#'
#' @section Errors:
#' \pkg{numDeriv} is in `Suggests` and every row needs it, so the function
#' stops when it is not installed. A `theta` outside the penalty's open bounds
#' or missing a hyperparameter is rejected before any check runs.
#'
#' @examples
#' # Second differences over five coefficients: rank 4, one null direction.
#' pen <- quadratic_penalty(crossprod(diff(diag(5))))
#' res <- check_penalty(pen)
#' all(res$status == "OK")
#'
#' # A separable penalty has fewer rows: no matrix, so no quadratic checks.
#' nrow(check_penalty(lasso_penalty(n_coef = 4), verbose = FALSE))
#'
#' # The validator earns its keep on a penalty that is wrong. Register the
#' # broken method on a subclass: registering on the real class would mutate
#' # the generic for the rest of the session.
#' Broken <- S7::new_class("Broken", parent = QuadraticPenalty)
#' bad <- do.call(Broken, S7::props(quadratic_penalty(diag(3))))
#' S7::method(penalty_gradient, Broken) <- function(pen, beta, theta, ...) {
#'   1.05 * S7::method(penalty_gradient, QuadraticPenalty)(pen, beta, theta, ...)
#' }
#' failed <- check_penalty(bad, verbose = FALSE)
#' failed[failed$status != "OK", c("check", "max_error")]
#'
#' @seealso [penalty_value()] and [penalty_gradient()] for the quantities
#'   checked, [penalty_kinks()] for the set the grid avoids,
#'   [check_abs_smoother()] for the same service on a smoother,
#'   [linkfunctions7::check_link()] and [distributions7::check_distrib()] for
#'   the siblings this follows.
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
    # The seed is fixed so that two runs report the same worst error: a
    # validator whose grid moves reports a different number every time. The
    # caller's own state is put back on exit, so a call in the middle of a
    # simulation leaves that simulation unchanged, and a caller who had drawn
    # nothing is left with no .Random.seed rather than with this one.
    old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
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
