# The Weighted Sum and Its Pseudo-Inverse

Assembles \\S(\lambda) = \sum_k \lambda_k P_k\\, takes one
eigendecomposition of it, and returns the matrix, its pseudo-inverse and
its log pseudo-determinant together. The value, the hyperparameter
derivatives and
[`penalty_logpdet()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.md)
all need some of these three, and this is the only place the
decomposition happens.

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

A list of three: `S`, the assembled symmetric matrix of side
`pen@n_coef`; `Sp`, its Moore-Penrose pseudo-inverse over the leading
`pen@p_rank` eigendirections, of the same side; and `logpdet`, a single
number, the sum of the logarithms of those eigenvalues. `Sp` and
`logpdet` are `NaN` where the decomposition cannot resolve the smallest
of those eigenvalues.

## Details

The pseudo-inverse is taken over the `p_rank` largest eigenvalues, that
rank having been fixed at construction from the components rather than
counted here. Selecting by rank instead of by a tolerance is what keeps
the answer steady when the parameters differ by many orders of
magnitude, which is exactly when a count would lose directions the
penalty still spans.

Selecting by rank is not on its own enough to make the answer
resolvable. A symmetric eigendecomposition computes its eigenvalues with
an absolute accuracy of order \\d\\\epsilon\lVert S\rVert\_{2}\\, so
where the parameters differ by enough orders of magnitude the smallest
eigenvalue spanning the range falls below that resolution and carries no
significant digit. The ordering between it and the null directions is
then rounding, and the sum of logarithms is `NaN` where the selected
value is negative and a plausible wrong number where it is positive.
Both are rejected: `Sp` and `logpdet` are `NaN` there. `S` is assembled
without a decomposition and is returned as it stands, which is why
[`penalty_hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_gradient.md)
stays finite.

## See also

[`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md),
[`penalty_logpdet.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
