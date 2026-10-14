# Quadratic Form of an Additive Penalty

Returns `TRUE`. A sum of quadratic forms is a quadratic form, so this
class has methods for the matrix, the rank, the null basis and the log
pseudo-determinant that
[`is_quadratic()`](https://statmodels7.github.io/penalties7/reference/is_quadratic.md)
gates, and a marginal criterion can read them. The matrix moves with one
hyperparameter per component, whereas the plain quadratic branch has one
scale, so the log pseudo-determinant is not linear in the logarithm of
any one parameter, and
[`check_penalty()`](https://statmodels7.github.io/penalties7/reference/check_penalty.md)
compares its gradient with `numDeriv` where the plain quadratic branch
checks a slope.

## Arguments

- pen:

  An
  [`AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/AdditivePenalty.md)
  object.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

`TRUE`.

## See also

[`penalty_matrix.AdditivePenalty()`](https://statmodels7.github.io/penalties7/reference/penalty_matrix.AdditivePenalty.md)
for the quantities that this predicate gates.
