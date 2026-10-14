# Drawing From a Separable Penalty

Draws from the parent family, one draw per coordinate for a univariate
parent and one per block for a multivariate parent, carried back through
the map.

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

A numeric vector of length `pen@n_coef`, `NA` where a diagonal map has a
zero entry and `NA` throughout under a map that is not diagonal.

## Details

The penalty is the parent's negative log-density applied to each
coordinate (or block) of \\D\beta\\, so a draw is
[`distributions7::distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.html)
at the hyperparameters, which are the parent's own parameters. Under a
multivariate parent the draw has one row per block and is flattened in
the order in which
[`dp_arg()`](https://statmodels7.github.io/penalties7/reference/dp_arg.md)
reads it.

A diagonal map is inverted by dividing, since the prior is on
\\d_j\beta_j\\. A coordinate whose \\d_j\\ is zero is not determined by
the prior and is `NA`, and the other coordinates are drawn as usual. A
square map whose off-diagonal entries are all zero is stored as a
diagonal map by
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md),
whatever its class, and is treated in the same way. Under any other map
every coordinate is `NA`.

## See also

[`penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.md),
[`distrib_penalty()`](https://statmodels7.github.io/penalties7/reference/distrib_penalty.md)
