# Drawing From a Separable Penalty

One draw per coordinate from the parent family, carried back through the
map.

## Arguments

- pen:

  A
  [`DistribPenalty()`](https://statmodels7.github.io/penalties7/reference/DistribPenalty.md)
  object.

- theta:

  Its hyperparameters, which are the parent's parameters.

- ...:

  Ignored.

## Value

A numeric vector of length `pen@n_coef`, all `NA` under a non-diagonal
map.

## Details

The penalty is the parent's negative log-density applied to each
coordinate of \\D\beta\\, so a draw is
[`distributions7::distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.html)
at the hyperparameters, which are the parent's own parameters. Under a
multivariate parent the draw is one row per block, flattened the way
[`dp_arg()`](https://statmodels7.github.io/penalties7/reference/dp_arg.md)
reads it.

A diagonal map is inverted by dividing, the prior being on
\\d_j\beta_j\\; any other map comes back `NA`, the coefficients not
being determined by a prior on a lower-dimensional image of them.

## See also

[`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md),
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
