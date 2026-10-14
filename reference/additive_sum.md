# The Weighted Sum and Its Log Pseudo-Determinant

Assembles \\S(\lambda) = \sum_k \lambda_k P_k\\ and returns it together
with its log pseudo-determinant and the two derivative blocks
\\\operatorname{tr}(S^{+}P_k)\\ and
\\\operatorname{tr}(S^{+}P_kS^{+}P_l)\\. The value, the hyperparameter
derivatives and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
of an
[`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
read these quantities from this function, which calls
[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md)
for the last three.

## Usage

``` r
additive_sum(pen, theta)
```

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- theta:

  The aligned hyperparameter list, as
  [`align_ptheta()`](https://statmodels7.github.io/penalties7/reference/align_ptheta.md)
  returns it.

## Value

A list of four: `S`, the assembled symmetric matrix of side
`pen@n_coef`; `logpdet`, a single number; `dlog`, a numeric vector of
\\\operatorname{tr}(S^{+}P_k)\\ in `pen@params` order; and `d2log`, the
square matrix of \\\operatorname{tr}(S^{+}P_kS^{+}P_l)\\. The last three
are `NaN` at a setting the transformation cannot resolve.

## A single decomposition of the sum

A single eigendecomposition of \\S(\lambda)\\, keeping its `p_rank`
largest eigenvalues, fails once the parameters spread apart. A symmetric
eigendecomposition computes its eigenvalues with an absolute accuracy of
order \\d\\\epsilon\lVert S\rVert\_{2}\\, so the smallest eigenvalue on
the range can fall below that resolution and carry no significant digit.
The sum of logarithms is then `NaN` if that eigenvalue is computed as
negative and an inaccurate number if it is computed as positive, and
neither case raises a warning.

## The similarity transformation

[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md)
partitions the components by the size of their contributions,
\\\lambda_k\max\_{ij}\lvert(P_k)\_{ij}\rvert\\, and rotates onto the
eigenvectors of the dominant group, the components whose contribution is
at least
[`additive_tol()`](https://statmodels7.github.io/penalties7/reference/additive_tol.md)
times the largest. With \\U\_{+}\\ spanning the range of that group and
\\U\_{0}\\ its kernel, \\A = U\_{+}^{\top}SU\_{+}\\, \\C =
U\_{+}^{\top}SU\_{0}\\, \\F = A^{-1}C\\ and \\M = U\_{0}^{\top}SU\_{0} -
C^{\top}F\\,

\$\$\log\mathrm{pdet}\\S = \log\lvert A\rvert + \log\mathrm{pdet}\\M,
\qquad S^{+} = U\_{+}A^{-1}U\_{+}^{\top} + GM^{+}G^{\top}, \qquad G =
U\_{0} - U\_{+}F.\$\$

Both identities are exact, because the null space of \\S\\ lies in the
span of \\U\_{0}\\. In the first, \\A\\ is the whole sum restricted to
the range of the dominant group, where the dominant group sets the
scale, and \\M\\ is a Schur complement of the size of the subordinate
components, so the decomposition of \\A\\ does not have to resolve the
spread between the two. The second writes every trace against \\S^{+}\\
as a trace against \\A^{-1}\\ plus a trace of the same kind against
\\M^{+}\\, so the derivatives are taken in the transformed coordinates
and never from a materialized \\S^{+}\\. \\M\\ is treated in the same
way at the next level. The recursion stops at a level where every
remaining component joins the dominant group, or where the rank of the
dominant group is zero or reaches the rank left, and that level takes
one eigendecomposition.

Two quantities are exactly zero and are dropped instead of computed. The
dominant group vanishes on \\U\_{0}\\, that subspace being its kernel,
so a computed value there would consist of rounding errors of order
\\\lambda\_{\mathrm{dom}}\epsilon\\, which would swamp the subordinate
terms; the blocks of each dominant component on \\U\_{0}\\ vanish for
the same reason. The part of \\M\\ that the next level does not carry as
components is likewise built from the reductions of the dominant
components and not obtained by subtracting the subordinate reductions
from \\M\\, since that subtraction would be a difference of two
quantities of the dominant size.

## The rank identity

\\\operatorname{tr}(S^{+}S) = r\\ exactly, so
\\\sum_k\lambda_k\operatorname{tr}(S^{+}P_k)\\ must equal the rank. The
identity needs no reference and is checked at every call. Where the
value returned by
[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md)
is not finite, where the identity is violated by more than
[`additive_check_tol()`](https://statmodels7.github.io/penalties7/reference/additive_check_tol.md)
times the rank, or where it cannot be evaluated because a trace
overflows, the point is treated as not resolvable in double precision,
and `logpdet`, `dlog` and `d2log` are `NaN`, without an error. `S` is
assembled without any decomposition and is returned as it is.
[`penalty_matrix()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md),
[`penalty_gradient()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
and
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
assemble the sum themselves and do not call this function, so they stay
finite at such a point unless the weighted sum itself overflows.

## See also

[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md)
for the transformation,
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
[`penalty_logpdet.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
