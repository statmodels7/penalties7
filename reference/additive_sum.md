# The Weighted Sum and the Determinant It Carries

Assembles \\S(\lambda) = \sum_k \lambda_k P_k\\ and returns it together
with its log pseudo-determinant and the two derivative blocks
\\\operatorname{tr}(S^{+}P_k)\\ and
\\\operatorname{tr}(S^{+}P_kS^{+}P_l)\\. The value, the hyperparameter
derivatives and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
all read some of these, and this is the only place the decomposition
happens.

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

## Why one decomposition is not enough

Taking a single eigendecomposition of \\S(\lambda)\\ and keeping its
`p_rank` largest eigenvalues fails once the parameters spread apart. A
symmetric eigendecomposition computes its eigenvalues with an absolute
accuracy of order \\d\\\epsilon\lVert S\rVert\_{2}\\, so the smallest
eigenvalue spanning the range falls below that resolution and carries no
significant digit: the sum of logarithms is then `NaN` where that value
comes back negative and a plausible wrong number where it comes back
positive. Measured on one anisotropic `te()` fit, 8419 calls of 12026
returned the second, with no warning, and the log pseudo-determinant was
out by as much as 618 against the exact asymptote.

## The similarity transformation

[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md)
partitions the components by the size they contribute, \\\lambda_k\lVert
P_k\rVert\\, and rotates onto the eigenvectors of the dominant group.
With \\U\_{+}\\ spanning that group's range and \\U\_{0}\\ its kernel,
\\A = U\_{+}^{\top}SU\_{+}\\, \\C = U\_{+}^{\top}SU\_{0}\\, \\F =
A^{-1}C\\ and \\M = U\_{0}^{\top}SU\_{0} - C^{\top}F\\,

\$\$\log\mathrm{pdet}\\S = \log\lvert A\rvert + \log\mathrm{pdet}\\M,
\qquad S^{+} = U\_{+}A^{-1}U\_{+}^{\top} + GM^{+}G^{\top}, \qquad G =
U\_{0} - U\_{+}F.\$\$

Both are exact. The first says no decomposition ever has to resolve the
spread, since \\A\\ carries the dominant group alone and \\M\\ what is
left; the second says every trace is a well-scaled trace against
\\A^{-1}\\ plus the same question one level down, so the derivatives are
taken in the transformed coordinates and never from a materialized
\\S^{+}\\. Measured, taking them in the original coordinates instead
loses one digit per order of magnitude of spread:
\\\lambda_1\operatorname{tr}(S^{+}P_1)\\ read 99537 where it is 15.

Two quantities are **exactly zero and are dropped rather than
computed**, and the accuracy rests on both. The dominant group vanishes
on \\U\_{0}\\, that subspace being its kernel, so computing it there
costs \\O(\lambda\_{\mathrm{dom}}\epsilon)\\ and swamps the subordinate
terms; and a dominant component's own \\U\_{0}\\ blocks vanish for the
same reason. What the next level does not carry is likewise built from
the dominant reductions rather than obtained by subtracting the
subordinate ones from \\M\\, which would be a difference of two
quantities of the dominant size whose difference is of the subordinate
one: measured, that spelling left a third component's log-scale gradient
at `5.0e-04` where it is exactly 8.

## What is rejected

\\\operatorname{tr}(S^{+}S) = r\\ exactly, so
\\\sum_k\lambda_k\operatorname{tr}(S^{+}P_k)\\ must be the rank. The
identity costs nothing, needs no reference and is checked at every call;
where it is violated by more than
[`additive_check_tol()`](https://statmodels7.github.io/penalties7/reference/additive_check_tol.md)
the point is not resolvable in double precision and `logpdet`, `dlog`
and `d2log` are `NaN`. `S` is assembled without any decomposition and is
returned as it stands, which is why
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
stays finite there.

## See also

[`additive_stable()`](https://statmodels7.github.io/penalties7/reference/additive_stable.md)
for the transformation,
[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
[`penalty_logpdet.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
