# An Additive Penalty Is Quadratic in the Coefficients

Answers `TRUE`. A sum of quadratic forms is a quadratic form, so this
branch has the matrix, the rank, the null basis and the log
pseudo-determinant that
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
gates, and a marginal criterion can read them. What is particular about
it is that the matrix moves with one hyperparameter per component, where
the plain quadratic branch has one scale, so the log pseudo-determinant
is not linear in any one of them and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
compares its gradient against `numDeriv` where the plain quadratic
branch reads a slope.

## Arguments

- pen:

  An
  [`additive_penalty()`](https://statmodels7.github.io/penalties7/reference/additive_penalty.md)
  object.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

`TRUE`.

## See also

[`penalty_matrix.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
for the quantities this admits.
