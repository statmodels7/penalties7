#' @include scad_mcp.R
NULL

#' @title Check a Penalty Numerically
#'
#' @description
#' Compares the derivatives that a penalty declares with numerical
#' references, and returns one row per comparison carrying the worst scaled
#' error and a pass or fail. The gradient is compared with a \pkg{numDeriv}
#' gradient of the value, and the Hessian with a numerical Jacobian of the
#' analytic gradient. Each hyperparameter gradient is compared with a
#' numerical derivative of the value in that hyperparameter, each entry of
#' the hyperparameter Hessian, diagonal or off-diagonal, with a numerical
#' derivative of the analytic hyperparameter gradient, and each mixed block
#' with a numerical Jacobian of the analytic coefficient gradient, so that no
#' reference is a difference of another difference. A quadratic penalty also has its
#' three-point identity, its log pseudo-determinant, its rank and its null
#' basis tested.
#' The function serves to verify the derivatives of a penalty written by the
#' user.
#'
#' @details
#' # The rows
#'
#' Two rows are produced for every penalty, three more for each
#' hyperparameter, one more for each pair of distinct hyperparameters, and
#' four or five more when [is_quadratic()] is `TRUE`:
#'
#' | row | compares |
#' |---|---|
#' | `gradient vs numDeriv` | [penalty_gradient()] against a numerical gradient of [penalty_value()] |
#' | `hessian vs numDeriv on the gradient` | [penalty_hessian()] against a numerical Jacobian of the analytic gradient, symmetrized |
#' | `grad_theta[p] vs numDeriv` | [penalty_grad_theta()] against a numerical derivative of the value in hyperparameter `p` |
#' | `hess_theta[p_p] vs numDeriv` | the diagonal entry of [penalty_hess_theta()] against a numerical derivative of the analytic `grad_theta` |
#' | `hess_theta[p_r] vs numDeriv` | the off-diagonal entry of [penalty_hess_theta()] for hyperparameters `p` and `r`, with `p` before `r` in `pen@params`, against a numerical derivative in `r` of the analytic `grad_theta` of `p` |
#' | `cross[p] vs Richardson on the gradient` | [penalty_cross()] against a numerical Jacobian of the analytic coefficient gradient in `p` |
#' | `quadratic three-point identity` | \eqn{\rho(2\beta) - 4\rho(\beta) + 3\rho(0)}, which vanishes when the value is a quadratic form with no linear term |
#' | `logpdet linear in log lambda with slope r` | raising `lambda` by a factor of \eqn{e} raises the value of [penalty_logpdet()] by the rank |
#' | `logpdet gradient vs numDeriv` | the `grad` element of `penalty_logpdet()` against a numerical derivative of its `value` |
#' | `rank vs eigenvalues of the matrix` | [penalty_rank()] against the number of eigenvalues of \eqn{S(\theta)}, or of the sum of the normalized components of an additive penalty, that exceed a relative tolerance times the largest |
#' | `null basis width vs n_coef minus rank` | the number of columns of [penalty_null_basis()] against `pen@n_coef` minus [penalty_rank()] |
#' | `null basis annihilates the matrix` | \eqn{S(\theta)N}, where \eqn{S(\theta)} is [penalty_matrix()] and \eqn{N} is [penalty_null_basis()] |
#'
#' The two log pseudo-determinant rows are alternatives, as described below.
#' The row of the annihilation appears only when the null basis has
#' columns. With
#' \eqn{p} hyperparameters the count is therefore \eqn{2 + 3p + p(p-1)/2} for
#' a penalty that is not quadratic, and four or five more for a quadratic
#' one. A lasso gives 5 rows, SCAD 9, a full-rank ridge 9, a quadratic
#' penalty over second differences 10, and a structured penalty over a
#' `3 x 3` log-Cholesky precision 39.
#'
#' # Which log pseudo-determinant row a quadratic branch gets
#'
#' The quadratic rows are produced for every branch for which
#' [is_quadratic()] returns `TRUE`: the three-point identity, one of the two
#' log pseudo-determinant rows, the row of the rank, the row of the width of
#' the null basis, and the row of the annihilation when the null basis has
#' columns. When `"lambda"` is among `pen@params`, the matrix is one scale
#' times a constant, the log pseudo-determinant is linear in
#' \eqn{\log\lambda} with slope the rank, and the row tests that identity.
#' Otherwise the matrix depends on several hyperparameters, as in a
#' structured or an additive penalty, and the row compares the reported
#' gradient with a numerical derivative of the value. A two-component
#' additive penalty whose sum has a null space therefore produces 14 rows:
#' the 9 rows of its two hyperparameters, the three-point identity, the
#' gradient of the log pseudo-determinant, the rank, the width of the null
#' basis and its annihilation. When the sum has full rank the row of the
#' annihilation is absent and the count is 13.
#'
#' The row of the rank counts the eigenvalues that exceed a relative
#' tolerance times the largest. For a quadratic or a structured penalty they
#' are the eigenvalues of \eqn{S(\theta)} at the `theta` of the check. For an
#' additive penalty they are the eigenvalues of the sum of the components,
#' each divided by its largest absolute entry, so that the row reproduces the
#' count of [additive_penalty()]: a count on \eqn{S(\theta)} drops the
#' smaller components once their scales differ by more than the tolerance,
#' as that function's page describes. The tolerance is the `tol` stored in a
#' quadratic or an additive penalty, the one given to its constructor, and
#' `1e-10` for a structured penalty, which stores none. The row's
#' `max_error` is the absolute difference between that count and
#' [penalty_rank()]. A stored rank that differs from the
#' count by one or more therefore fails the row of the rank, and a null
#' basis whose number of columns differs from `pen@n_coef` minus the stored
#' rank fails the row of the width.
#'
#' # What is not checked
#'
#' The check covers the value and its derivatives on the parameter scale, and
#' the quantities that a marginal criterion reads. It does not test the
#' `hess` element of [penalty_logpdet()], the `scale = "link"` option of the
#' hyperparameter derivatives, [is_proper()], [penalty_prox()],
#' [penalty_prox_spec()], [penalty_dhessian()], [penalty_d2hessian()],
#' [penalty_dcross()], [penalty_dhessian_beta()], [penalty_d2hessian_beta()],
#' [penalty_dhessian_beta_theta()], [penalty_draw()] or [penalty_readable()].
#' It calls [penalty_kinks()] to place the coefficients but does not test the
#' kinks themselves.
#'
#' # Where the coefficients are placed
#'
#' A penalty with a kink has no derivative there, so a numerical reference
#' that straddles one measures the kink and not the formula. With `beta` left
#' at `NULL`, `0.033` is added to every coefficient, up to fifty times, until
#' every coordinate of \eqn{D\beta} lies more than `0.05` from every kink that
#' [penalty_kinks()] returns. A draw that already clears the kinks is used
#' unchanged, and after fifty steps the last draw is used whether or not it
#' clears them.
#'
#' @param pen A [penalty()] object, of any branch.
#' @param beta A numeric coefficient vector of length `pen@n_coef`. `NULL`,
#'   the default, draws one from `rnorm(sd = 1.3)` rounded to two places and
#'   shifted by `0.11`, then moves it away from the kinks as above.
#'   The draw is taken from a fixed seed, so the report is the same on two
#'   runs, and the caller's own `.Random.seed` is restored on exit (or removed,
#'   if there was none), so a call in the middle of a simulation leaves that
#'   simulation unchanged.
#' @param theta A named list of hyperparameter values, or a named numeric
#'   vector carrying the same. `NULL`, the default, places each hyperparameter
#'   six tenths of the way across its own bounds, reading an infinite lower
#'   bound as `-1` and an infinite upper bound as two above the lower. That
#'   gives `lambda = 1.2` for a hyperparameter bounded on
#'   \eqn{(0, \infty)}, `alpha = 0.6` for the elastic net, `a = 3.2` for SCAD
#'   and `gamma = 2.2` for MCP.
#' @param tol The scaled error at or above which a row is reported as
#'   `FAILED`. A single positive number, `1e-6` by default. A row whose error
#'   is not finite is also reported as `FAILED`, and so is a row whose
#'   numerical reference cannot be computed.
#' @param verbose `TRUE`, the default, prints the table without row names.
#'   The result is returned invisibly either way.
#'
#' @return A data frame with one row per check and three columns: `check`
#'   (character, the row's name as tabulated above), `max_error` (numeric) and
#'   `status` (character, `"OK"` or `"FAILED"`). Returned invisibly. For the
#'   rows that compare with a numerical reference, `max_error` is the worst
#'   absolute difference divided by `max(1, max(abs(reference)))`. For the
#'   three-point identity it is the absolute value of the combination divided
#'   by `max(1, abs(rho(beta)))`; for the slope of the log pseudo-determinant,
#'   the absolute difference between the increase and the rank; for the rank
#'   and for the width of the null basis, the absolute difference between the
#'   two counts; for the annihilation by the null basis, the largest absolute
#'   entry of \eqn{S(\theta)N} divided by `max(1, max(abs(S)))`. A row whose
#'   reference signals an error, a numerical derivative or the eigenvalues
#'   of the rank row, has `max_error` equal to `NaN`.
#'
#' @section Errors:
#' \pkg{numDeriv} is in `Suggests` and the derivative rows need it, so the
#' function signals an error when it is not installed. A `theta` outside the
#' penalty's open bounds or missing a hyperparameter is rejected before any
#' check runs. The analytic methods are evaluated at `beta` and `theta`
#' outside the numerical references, so an error in [penalty_value()],
#' [penalty_gradient()], [penalty_hessian()], the hyperparameter derivatives
#' or the quadratic quantities is signaled as it arises. An error inside a
#' numerical reference is not signaled: \pkg{numDeriv} stops when the
#' function it differentiates returns a value that is not finite, and the
#' row whose reference stopped is then reported with `max_error` equal to
#' `NaN` and status `FAILED`, while the other rows are computed as usual.
#'
#' @examples
#' # First differences over five coefficients: rank 4, one null direction.
#' pen <- quadratic_penalty(crossprod(diff(diag(5))))
#' res <- check_penalty(pen)
#' all(res$status == "OK")
#'
#' # A separable penalty has fewer rows: no matrix, so no quadratic checks.
#' nrow(check_penalty(lasso_penalty(n_coef = 4), verbose = FALSE))
#'
#' # A penalty with a wrong gradient fails the rows that read the gradient.
#' # The broken method is registered on a subclass, because registering it on
#' # the real class would change the generic for the rest of the session.
#' Broken <- S7::new_class("Broken", parent = QuadraticPenalty)
#' bad <- do.call(Broken, S7::props(quadratic_penalty(diag(3))))
#' S7::method(penalty_gradient, Broken) <- function(pen, beta, theta, ...) {
#'   1.05 * S7::method(penalty_gradient, QuadraticPenalty)(pen, beta, theta, ...)
#' }
#' failed <- check_penalty(bad, verbose = FALSE)
#' failed[failed$status != "OK", c("check", "max_error")]
#'
#' @seealso [penalty_value()] and [penalty_gradient()] for the quantities
#'   checked, [penalty_kinks()] for the set that the coefficients avoid,
#'   `numericals7::check_abs_smoother()` for the same check on a smoother,
#'   [linkfunctions7::check_link()] and [distributions7::check_distrib()] for
#'   the corresponding checks of a link and of a distribution.
#' @export
check_penalty <- function(pen, beta = NULL, theta = NULL, tol = 1e-6,
                          verbose = TRUE) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("check_penalty() requires the numDeriv package.", call. = FALSE)
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
  # a numerical reference that stops (numDeriv does on a non-finite value)
  # gives NaN, which add() reports as FAILED; the analytic methods are called
  # outside it at beta and theta, so their own errors are still signaled
  reference <- function(expr) tryCatch(expr, error = function(e) NaN)

  penalty_value(pen, beta, theta)
  g <- penalty_gradient(pen, beta, theta)
  gn <- reference(numDeriv::grad(function(b) penalty_value(pen, b, theta), beta))
  add("gradient vs numDeriv", rel(g, gn))

  H <- penalty_hessian(pen, beta, theta)
  Hn <- reference(numDeriv::jacobian(function(b) penalty_gradient(pen, b, theta),
                                     beta))
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
        rel(gt[[p]], reference(numDeriv::grad(fn, theta[[p]]))))
    gfn <- function(v) {
      th <- theta; th[[p]] <- v
      penalty_grad_theta(pen, beta, th)[[p]]
    }
    add(paste0("hess_theta[", p, "_", p, "] vs numDeriv"),
        rel(ht[[paste0(p, "_", p)]],
            reference(numDeriv::grad(gfn, theta[[p]]))))
    cfn <- function(v) {
      th <- theta; th[[p]] <- v
      penalty_gradient(pen, beta, th)
    }
    add(paste0("cross[", p, "] vs Richardson on the gradient"),
        rel(cr[[p]],
            reference(as.numeric(numDeriv::jacobian(cfn, theta[[p]])))))
  }
  # the mixed entries of the hyperparameter Hessian: the analytic gradient in
  # the first hyperparameter of the pair, differentiated in the second
  prs <- ptheta_pairs(pen@params)
  for (key in names(prs)) {
    pr <- prs[[key]]
    if (pr[1] == pr[2]) next
    mfn <- function(v) {
      th <- theta; th[[pr[2]]] <- v
      penalty_grad_theta(pen, beta, th)[[pr[1]]]
    }
    add(paste0("hess_theta[", key, "] vs numDeriv"),
        rel(ht[[key]], reference(numDeriv::grad(mfn, theta[[pr[2]]]))))
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
        ref <- reference(numDeriv::grad(function(v) {
          th <- theta; th[[p2]] <- v
          penalty_logpdet(pen, th)$value
        }, theta[[p2]]))
        worst <- max(worst, abs(lp$grad[[p2]] - ref) / max(1, abs(ref)))
      }
      add("logpdet gradient vs numDeriv", worst)
    }
    M <- penalty_matrix(pen, theta)
    nb <- penalty_null_basis(pen)
    # the stored rank against the eigenvalues of the matrix, counted by the
    # relative rule of the constructors at the tolerance the penalty stores
    # (the default one for a branch that stores none), and the width of the
    # stored null basis against the stored rank. An additive penalty is
    # counted on its normalized components, as additive_penalty() counts it:
    # S(theta) loses the smaller components once their scales differ by more
    # than the tolerance.
    r_tol <- if (S7::S7_inherits(pen, QuadraticPenalty) ||
                 S7::S7_inherits(pen, AdditivePenalty)) pen@tol else 1e-10
    Md <- if (S7::S7_inherits(pen, AdditivePenalty)) {
      as.matrix(Reduce(`+`, lapply(pen@mats, function(P) P / max(abs(P)))))
    } else {
      as.matrix(M)
    }
    count <- reference({
      ev <- eigen((Md + t(Md)) / 2, symmetric = TRUE, only.values = TRUE)$values
      sum(ev > r_tol * max(ev, 0))
    })
    r <- penalty_rank(pen)
    add("rank vs eigenvalues of the matrix", abs(count - r))
    add("null basis width vs n_coef minus rank",
        abs(ncol(nb) - (q - r)))
    if (ncol(nb)) {
      add("null basis annihilates the matrix",
          max(abs(M %*% nb)) / max(1, max(abs(M))))
    }
  }

  out <- do.call(rbind, checks)
  if (verbose) print(out, row.names = FALSE)
  invisible(out)
}
