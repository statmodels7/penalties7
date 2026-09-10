# penalties7 0.22.0

* `beta_quadratic()` reads a MULTIVARIATE parent's response Hessian between
  readings rather than between the entries of one of them. The predicate
  decides whether the penalized objective is quadratic in the coefficients,
  hence whether the inner problem is one linear solve, and it asked whether
  every entry of `distrib_hess_y()` equals the first. For a parent over more
  than one coordinate that compares an off-diagonal with a diagonal:
  measured on the multivariate gaussian a covariance class declares, whose
  Hessian is the constant -Sigma^-1 reading -1, 0, 0, -1, the predicate
  answered FALSE where the prior is exactly quadratic, at every dimension
  from two to four. It answers TRUE now, and a multivariate Student t, whose
  Hessian carries the (nu + p)/(nu + q) reweighting and so does move with
  the response, still answers FALSE -- which is what says the change is not
  simply answering TRUE for anything multivariate.

* The probe is four POINTS of `pen@block` coordinates each rather than four
  VALUES. A p-variate parent reads a point as p numbers, so four values laid
  out with p columns are two points at p = 2 or 3 and ONE from p = 4 up, and
  a single reading takes the branch that answers without comparing anything.
  Measured on a multivariate Student t prior, where the truth is FALSE at
  every width, four values answer FALSE at p = 2 and 3 and TRUE at 4 and 5;
  four points answer FALSE at all of them. Four values also recycled at an
  odd width, so `beta_quadratic()` of a three-variate parent emitted a
  `data length [4] is not a sub-multiple` warning from inside a fit.

* A UNIVARIATE parent is untouched, and the offset that spaces the points is
  exactly zero at a block width of one, so such a parent is probed at the
  same four values it always was. Measured against the form this release
  replaces, a Gaussian prior, a Student t prior and a Laplace prior answer
  TRUE, FALSE and TRUE under both.

* `test-marginal-derivatives.R` had no multivariate case at all, which is
  how this survived. It carries three now, and each of the two halves has
  its own negative control written out as the form it replaces: restoring
  the entrywise comparison fails the gaussian block three times, once per
  width, and restoring the four-value probe fails the Student t block and
  emits the recycling warning with it. The univariate block stays green
  under both, which is what says the change is confined to the multivariate
  branch.

# penalties7 0.21.0

* `penalty_draw()` draws from a penalty read as a prior. A penalty is a
  negative log-density -- which is why the normalizing constant is kept --
  so it is something one can draw from, and a simulation that carries a
  random effect, a ridge or a lasso should take those coefficients from
  that distribution rather than from a normal chosen by whoever wrote the
  simulation.

* A separable penalty draws coordinatewise from its own parent family, so a
  Gaussian prior gives Gaussian effects, a Laplace prior gives Laplace ones
  and a Student t prior heavy-tailed ones. Checked against the parent's own
  distribution function, which the draw does not read: a Laplace draw passes
  that test and fails the Gaussian of the same variance at 1e-6.

* A quadratic, additive or structured penalty is a Gaussian prior with
  precision S(theta), and two of its shapes are cheap: a diagonal S, which
  is a ridge and a Demmler-Reinsch smooth, is drawn coordinatewise at
  1/sqrt(S_jj), and a full-rank S through one Cholesky factor. Measured
  against the inverse of the precision itself, the sampled covariance agrees
  to 0.03 over 60000 draws.

* A coordinate the prior does not determine comes back `NA` rather than
  zero: the unpenalized direction of a smooth, every coordinate of SCAD and
  MCP, which are improper by construction, a deficient non-diagonal
  precision, and a separable penalty under a general map, where the prior
  leaves the coefficients underdetermined. Zero would be indistinguishable
  from a prior that concentrates there, and a caller filling the rest by a
  rule of its own could not tell which coordinates those are.

# penalties7 0.20.0

* `structured_penalty()` reads its structure as the prior's **precision**,
  always. It no longer reads `parameters7`'s `role`, which that package
  removes at 0.19.0, and it no longer rejects a structure that declares
  itself usable as either -- there is nothing left to declare.

  The prior is a negative log-density and the matrix it needs is the one in
  the quadratic form, so the side is a property of the construction rather
  than a choice offered to the caller. The book has said so in
  `sec-penalty-structured` since it was written.

* A caller who wants the structure on the covariance side writes
  `structured_penalty(parameters7::inverse_of(s))`. Nothing becomes
  inexpressible: `inverse_of()` carries the chain rule for an inverse to
  fourth order inside the structure, where the log-determinant stays exact,
  and for `ar1()` and `autoregressive()` the inverse has a name of its own.
  The twin that made the sign verifiable -- the covariance reading at Sigma
  against the precision reading at Sigma^-1, agreeing to the last bit -- is
  kept, rewritten through the wrapper.

* `struct_is_cov()` is gone and the four helpers `struct_omega()`,
  `struct_d1()`, `struct_d2()` and `struct_logdet()` no longer transport: a
  structure meant for the other side is a different structure and supplies
  its own arrays. `penalty_name` drops the role it used to carry, reading
  `structured [ar1]` rather than `structured [ar1, precision]`.

* A rank-deficient structure is admitted without qualification, giving the
  improper prior the log pseudo-determinant is written for. The refusal that
  used to sit here for the covariance reading now sits in `inverse_of()`,
  where it belongs: a singular matrix has no inverse.

# penalties7 0.19.0

* `is_quadratic()` answers `TRUE` for an `additive_penalty()`. It inherited
  `FALSE` from the base class while `penalty_matrix()`, `penalty_rank()`
  and `penalty_logpdet()` all answered for that branch and
  `penalty_null_basis()` rejected with a message pointing at a predicate
  that would have told the reader nothing. A sum of quadratic forms is a
  quadratic form.

* `penalty_null_basis()` has a method there. The null space of the sum is
  the intersection of the components', so it does not move with the
  hyperparameters, and the constructor already decomposed the normalized
  components to read the rank; what changed is that it keeps the vectors.

* `penalty_logpdet()` on that branch returns its gradient and Hessian as
  named lists keyed by hyperparameter and by pair, which is the shape the
  quadratic and structured branches use. It returned an unnamed numeric
  vector and a square matrix. Nothing outside this package read either:
  every consumer routed on `is_quadratic()` first and never arrived.

* `check_penalty()` therefore runs its three quadratic rows on an additive
  penalty, which it did not. A two-component one goes from 8 rows to 11,
  gaining the three-point identity, the log pseudo-determinant's gradient
  against `numDeriv` and the null basis against the matrix, and its rank,
  matrix and log pseudo-determinant stop being untested.

* The proximal route is unchanged. `has_prox()` asks whether
  `penalty_prox()` is registered on the base class before it asks
  `is_quadratic()`, and this branch registers none, so it still answers
  `FALSE`.

# penalties7 0.18.0

* `check_penalty()` leaves the caller's random stream as it found it. With
  `beta` at its default the function called `set.seed(7)` and never restored
  what was there, so a call inside a simulation silently changed the
  simulation. Measured, `set.seed(99); runif(1)` gave 0.5847119 and the same
  with a `check_penalty()` call in between gave 0.3400624; passing `beta` skipped
  the draw and left the stream alone, which is what identified the branch.

  The seed stays fixed, so a validator still reports the same worst error on two
  runs, and the draw itself is unchanged, so no reported number moves: the
  default report is identical to the one taken at the beta that seed produces.
  What is new is the restore, through `on.exit()`, so it happens even when a
  check signals. A caller who had drawn nothing is left with no `.Random.seed`
  rather than with this one.

  This is the same habit `parameters7::check_parameter()` had, fixed in
  parameters7 0.13.0.

# penalties7 0.17.0

* A penalty with no free hyperparameters answers on the whole public surface. A
  fully known prior is a legitimate object: `distributions7::fixed()` documents
  `n_params = 0` as legal, and it is what a caller builds to hold a penalty at a
  value the outer search must not touch. Eight generics answered for one and
  three stopped, for two unrelated reasons.

  `ptheta_pairs()` built its names with `paste0(params, "_", params)`, which
  recycles a zero-length argument against the length-one literal and gives the
  single string `"_"` while the list of pairs is empty, so `setNames()` raised
  `'names' attribute [1] must be the same length as the vector [0]` three frames
  below the call, naming neither the penalty nor the argument. That stopped
  `penalty_hess_theta()` and, through it, `check_penalty()`, so such a penalty
  could not be validated at all. It returns the empty named list now, which is
  the shape `penalty_grad_theta()` and `penalty_cross()` already had for the
  same case.

  `penalty_prox()` read its hyperparameter as `theta[[which(pen@params ==
  "sigma")]]`, and `which()` of an empty comparison is `integer(0)`, so the
  subscript raised *attempt to select less than one element in get1index*. A
  held parameter is not absent, only kept elsewhere: `.prox_param()` takes it
  from `theta` when the penalty carries it free and from the parent's
  `fixed_params` when `fixed()` holds it. Measured, the three closed-form
  families now return exactly what the same penalty with the parameter free
  returns at that value, to the bit, for the operator and for the piecewise
  table the compiled coordinate descent reads; `penalty_value()` already agreed.
  A parent carrying no such parameter is named in the message.

  This also removes a contradiction: `has_prox()` answered `TRUE` for a penalty
  whose `penalty_prox()` stopped.

# penalties7 0.16.0

* New class `abs_smoother`: a smooth replacement `s(u)` for `|u|` carrying
  its derivatives in `u` up to order five as FUNCTIONS, so a piecewise
  smoother's branches are ordinary code rather than an expression
  `stats::deriv` cannot read. The contract is only `s`: the smooth sign is
  `s'`, the smooth step `(1 + s')/2` and the smooth hinge `(u + s)/2`
  follow by composition, which is what lets one object serve the smoothed
  break-point terms now and a smoothed kinked penalty later.

* Three instances: `smooth_probit()` (the recommended default: gaussian
  tails, and the exact convolution identity `tau_true^2 = tau^2 - h^2`
  declared as its `tau_correction`), `smooth_hyperbolic()`
  (`sqrt(u^2 + c)`, polynomial tails, no correction) and
  `smooth_quintic()` (exact outside `[-h, h]`, `C^3` at the seam).

* `smoother_deriv()`, `smoother_width()` (a `NULL` width is resolved by
  the consumer at build, from the covariate's spacing) and
  `smoother_width_floor()`, whose floor is derived from the expression
  that binds -- the Jacobian column carries `s''(0)/2 ~ 1/h` against
  columns of order the range `D`, so holding the design's condition below
  `eps^-1/2` gives `h >= sqrt(eps) * D` -- rather than chosen.

* `check_abs_smoother()`, the sibling of `check_penalty()` for
  user-written smoothers: the structural properties (even, odd bounded
  first derivative, convex, matching `|u|` in the tails) and each order
  against one numerical differentiation of the analytic order below it.

# penalties7 0.15.0

* `quadratic_penalty()` takes `blocks`: the penalty of `I_m (x) P`, built
  from `P` alone. What the constructor needs from its matrix is the rank, the
  log pseudo-determinant and a basis of the null space, and all three follow
  from one block -- the eigenvalues of `I_m (x) P` are `P`'s repeated `m`
  times -- so the assembled matrix is never formed or decomposed.

* Measured at `m = 200` over a basis of ten: 5.65 s to build from the
  assembled matrix against 0.010 s from the block, 565 times, of which 4.50 s
  was the eigendecomposition alone. The stored matrix is sparse besides, 0.03
  MB against 32.00, which follows rather than being the point. Rank, log
  pseudo-determinant, value, gradient, Hessian, both theta derivatives, the
  mixed block and the null space all agree with the penalty built from the
  assembled matrix, the gradient and the Hessian to exactly zero.

* It does not combine with `map`: a map mixes the blocks, and `D'(I (x) P)D`
  is block diagonal with a DIFFERENT block each, which is not the structure
  `blocks` names.

* This is the one branch whose `penalty_hessian()` is not a base matrix, and
  deliberately: it exists to avoid the assembled form, so returning that form
  would defeat it. `penalty_dhessian()` follows, `unclass()` having done
  nothing to an S4 object. A consumer coerces where the two kinds meet.

* The same identity is what `parameters7::kron_identity()` uses on the other
  side of the toolkit, for the covariance of grouped random effects.

# penalties7 0.14.1

* The blockwise marginal pieces are exercised on a Student t parent as well,
  whose response Hessian depends on the observation: `dp_blockdiag()` already
  carried both shapes, and the three quantities agree with one difference of
  the analytic quantity below each.

# penalties7 0.14.0

* `penalty_readable()` says what a penalty's hyperparameters are ABOUT, where
  they are coordinates of a chart rather than the quantities themselves: a
  correlated prior's are the logarithms of a Cholesky diagonal and the entries
  below it, and what the prior describes is the standard deviations and the
  correlations of the effects. It is the same distinction
  `parameters7::param_readable()` makes for a matrix parameter and
  `modelterms7::term_readable()` for a fitted term. The base method returns
  `NULL`, which says the hyperparameters ARE the quantities -- the honest
  answer for every other branch, a smoothing parameter, a rate and a shape
  each being read on their own scale already.

* `penalty_dhessian()`, `penalty_d2hessian()` and `penalty_dcross()` carry the
  blockwise reading, so a marginal criterion estimates a correlated prior's
  matrix exactly rather than refusing. Each is one difference of the analytic
  quantity below it against numDeriv, at p = 2 and 3.

# penalties7 0.13.0

* `distrib_penalty()` reads its parent one BLOCK at a time, the block being
  the parent's dimension. A univariate parent gives blocks of one, which is
  the separable penalty unchanged; a p-variate parent gives blocks of p,
  which is a prior letting the coordinates of one block depend on each other
  while the blocks stay independent -- the effects of one group of a random
  effect. The Hessian is block diagonal rather than diagonal, the mixed block
  comes from the parent's `distrib_cross_y()`, and nothing else moves.
  Against a hand-written negative log multivariate normal density the value
  agrees to 2e-15, and every derivative agrees with numDeriv.

* A blockwise parent has no proximal operator, and `has_prox()` says so
  before a fitting layer routes a block to a scheme that cannot solve it: the
  operator acts one coordinate at a time and a correlated block does not
  separate. `penalty_prox_spec()` returns NULL there for the same reason, and
  a kink is a point of a scalar argument, so `distrib_kinks()` of a
  multivariate parent is empty and says so rather than inheriting the
  univariate route by accident.

* `structured_penalty()` honors the structure's `role`. A structure declared
  `"covariance"` is read as the covariance of the prior, and the quantities
  below are the same arithmetic at the precision it implies, transported once
  by the chain rule for an inverse. A structure that declares `"either"` is
  now REJECTED: the two readings differ in the sign of the log-determinant
  term, so a guess would give a fit converging to a different matrix without
  saying so, and the property has been carried, validated and never read
  since it was added. A covariance of deficient rank is rejected too, having
  no inverse; a precision of deficient rank is the improper prior the log
  pseudo-determinant is written for and is still admitted. The two readings
  are pinned against each other at Sigma and Sigma^-1, where they agree to
  the last bit.

# penalties7 0.12.0

* `ridge_penalty()` is written by its PRECISION and built on the quadratic
  branch, so its hyperparameter is `lambda` and a larger value shrinks
  harder -- as it already did for the lasso, whose `laplace2` carries a
  rate. The Gaussian prior written by its precision IS the quadratic
  penalty at the identity, the same value to the last bit (5.177374 against
  5.177374 on three coefficients), so what the move removes is a second
  name for one number rather than a construction: the separable twin is
  still there, and the tests pin the two against each other as before, read
  the other way round.

  Two consequences. A random effect's hyperparameter is a precision, so the
  variance component a reader wants is `1/sqrt(lambda)`. And
  `penalty_prox_spec()` has no entry for the quadratic branch, so the ridge
  no longer has a piecewise table; nothing loses a route, a smooth penalty
  never reaching the compiled coordinate descent, and `penalty_prox()` is a
  linear solve there.

# penalties7 0.11.0

* `penalty_prox_spec()` survives a diagonal map, where it returned `NULL`
  under any map at all.

  The table is the route a compiled coordinate descent takes, so a
  standardized penalty that lost it would have fallen back on the general
  proximal operator and paid an R call per coordinate per sweep -- the
  half of the seam 0.10.0 left open. The same change of variable carries
  the table across: reading it at `d v` with the step `t d^2` and dividing
  back divides the cuts and the intercepts by `|d|` and leaves the slopes
  alone, the slope multiplying a point that was scaled and then divided.
  The operator is odd, so only the magnitude of `d` enters.

  Pinned against `penalty_prox()` under the same map, coordinate by
  coordinate and across every breakpoint including the breakpoints
  themselves, on all five separable families: 1e-12. A map that mixes
  coordinates still returns `NULL`.

  The convexity condition of SCAD and MCP is now tested on the scaled step,
  so it is `t < (a - 1)/d^2` and `t < gamma/d^2`; a step admissible under
  the identity map and not under a map that stretches returns `NULL` rather
  than a table of a non-convex subproblem.

# penalties7 0.10.0

* A DIAGONAL map keeps the proximal operator, where any map used to lose it.

  A diagonal `D` rescales each coordinate on its own, so a separable penalty
  under one is still separable and the operator follows from the identity one
  by a change of variable: with `u = d b`,

      argmin_b (b - v)^2/(2t) + rho(d b)
        = argmin_u (u - d v)^2/(2 t d^2) + rho(u),   b = u/d

  so it is the same closed form read at the scaled point with the step scaled
  by `d^2`, divided back. `has_prox()` answers TRUE there, and a map that is
  not diagonal is still rejected, mixing coordinates being the
  generalized-lasso problem rather than a different formula.

  This is what standardization is. Penalizing a column divided by its own
  spread is penalizing `rho(s_j beta_j)`, so the scaling never has to touch
  the design: the sparsity of a block survives it, and the centring that would
  destroy that sparsity is not needed at all where an intercept is free.
  Measured against a coordinate-by-coordinate minimization that shares no
  arithmetic with the closed forms, on the four separable penalties: 4.0e-06,
  4.7e-06, 5.0e-06 and 4.8e-06, which is the reference grid's own resolution.

  ⚠️ The convexity condition of SCAD and MCP binds on the SCALED step, so it
  becomes `t < (a - 1)/max(d^2)` and `t < gamma/max(d^2)`. The existing guards
  catch it unchanged -- with `d` up to 10 and `t = 0.3` the effective step is
  30 and is refused -- but a caller who standardizes takes shorter steps.

* `as_map()` keeps a map that is already a `Matrix`. The constructors called
  `as.matrix()` on it in five places, which densified a diagonal map into the
  `q x q` matrix it exists to avoid: that coercion, not the operator, was the
  real obstacle to carrying a scaling without giving up sparsity.

# penalties7 0.9.0

* `penalty_prox_spec()` describes the scalar proximal operator of a
  separable penalty as an odd piecewise linear table, so a compiled loop
  can apply it without knowing which family it came from. Every closed
  form the package carries has that shape: the soft threshold is two
  pieces, the elastic net two, MCP three, SCAD four, and a Gaussian prior
  one. The step is a vector, one per coefficient, because in a coordinate
  descent the step of coordinate j is 1/sum(w x_j^2) and does not move
  while the working weights are held.

  The operator is applied once per coordinate per sweep at a point that
  moves every time, so a compiled loop calling back for it would spend the
  gain on the calls. Passing the numbers instead keeps the mathematics in
  the penalty and leaves the kernel naming no family.

  `prox_apply()` evaluates a table in R, and the table is pinned against
  `penalty_prox()` itself across every breakpoint and at the breakpoints
  exactly, at three step lengths and five families. A penalty with no such
  description -- a quadratic under a general matrix, an operator that is a
  root rather than a formula, a parent not centered where the quadratic
  pull is, a step past the convex region of SCAD or MCP -- returns `NULL`.

# penalties7 0.8.0

* `distrib_penalty()` derives its kink set from the parent instead of
  defaulting to none. A penalty built by hand from a non-smooth family --
  `distrib_penalty(fixed(laplace_distrib(), mu = 0))`, which is the lasso --
  declared itself differentiable everywhere, so a model layer reading
  `penalty_kinks()` put it in the scheme for the opposite property. The
  shipped instances passed `kinks` explicitly and were never affected.

  `distrib_kinks()` takes the candidates from `params_smooth` crossed with
  what `fixed()` holds, a location that is not smooth being a kink in the
  argument at the value it is held at, and then measures each one by
  comparing the one-sided derivatives of the log-density across it.
  Inferring alone would put a kink on any family whose non-smooth parameter
  is not a location; measuring alone would need somewhere to look. Nothing
  is taken from a parameter that is free, its value being whatever the
  hyperparameters say at the time.

  Passing `kinks` still overrides, including `numeric(0)` to declare there
  are none.

# penalties7 0.7.1

* The hyperparameters may be given as a named numeric vector as well as
  the documented list, the alignment converting one to the other. The
  branches had split on how they read `theta`: `[[` accepts both shapes
  and `$` accepts only the list, so a caller passing a vector reached
  the quadratic and separable branches and stopped inside `scad()` and
  `mcp()` with "$ operator is invalid for atomic vectors", three frames
  down and naming neither the argument nor the penalty.

# penalties7 0.7.0

* `penalty_d2hessian()` and `penalty_dcross()` read the parent's
  `distrib_grad_y_hess()` and `distrib_hess_y_hess()` for the separable
  branch instead of differencing its first-order components. With a
  gaussian parent -- every ridge, every random effect -- the branch is
  now exact: measured, the second derivative of `I/sigma^2` comes back
  as `6I/sigma^4` to 1e-13 where the difference gave 1e-10.

# penalties7 0.6.0

* `penalty_dhessian()`, `penalty_d2hessian()` and `penalty_dcross()`
  are what a marginal criterion asks of a penalty: the hyperparameter
  derivatives of the coefficient Hessian and of the mixed block. Every
  branch answers -- the quadratic and additive ones from their
  components, the structured one from the matrix parameter's `param_d1`
  and `param_d2`, the separable one from the parent's response surface
  -- so a penalty is estimable by REML or ML whatever its shape, and one
  with a kink rejects by name. `beta_quadratic()` reports whether the
  third derivative in the coefficients is zero, asking the parent
  whether its response curvature depends on the response rather than
  recognizing a family by name.

# penalties7 0.5.0

* elasticnet_penalty(): the elastic net as the same construction as
  ridge and lasso, a separable penalty over a fixed() family -- here
  distributions7::enet_distrib(), the product of the Laplace and the
  Gaussian at zero, normalized. Its hyperparameters are the overall
  rate lambda and the mixing weight alpha, and the normalizing
  constant depends on both, which is what makes them estimable by a
  marginal criterion.
* penalty_prox() gains the elastic-net closed form, the soft
  threshold of the Laplace part followed by the shrinkage of the
  Gaussian one.

# penalties7 0.4.0

* additive_penalty(): a sum of quadratic penalties with a smoothing
  parameter on each, which is what an anisotropic tensor smooth needs.
  The rank is fixed at construction from the components stacked and
  normalized, since the null space of the sum is the intersection of
  theirs and does not move with the parameters -- a count taken from
  the assembled matrix does.

# penalties7 0.3.0

* penalty_prox(): the proximal operator, closed form for the quadratic
  and structured branches (one linear solve), for the Gaussian and
  Laplace instances, and for SCAD and MCP over their piecewise
  regions; any other separable penalty is solved coordinatewise from
  its response derivative. has_prox() asks before calling.

# penalties7 0.2.0

* The structured quadratic prior: `structured_penalty()` takes a
  parameters7 matrix_parameter as the PRECISION, so the hyperparameters
  reach every entry of the matrix. The free vector is unconstrained by
  construction, so every link is the identity -- the flattening convention
  of the multivariate families. Every derivative comes from the
  structure's own contract (param_d1/param_d2, the logdet derivatives),
  and the marginal pieces answer through it, so is_quadratic() is TRUE.
  At a zero free vector the log-Cholesky prior IS the plain ridge, pinned
  at machine precision. check_penalty() gains the structured logpdet
  check, its lambda-slope identity now gated on the branch that has a
  lambda.

# penalties7 0.1.0

* First release: penalties as S7 objects, rho(D beta; theta), with the
  value, the exact derivatives in the coefficients and the
  hyperparameters, the mixed block, and the declared kink set. Three
  branches: `quadratic_penalty()` (rank, null basis and log
  pseudo-determinant fixed at construction, the pieces a marginal
  criterion consumes), `distrib_penalty()` (a univariate distributions7
  log-density applied coordinatewise, with `ridge_penalty()`,
  `lasso_penalty()` and `heavy_penalty()` as named instances), and
  `scad_penalty()`/`mcp_penalty()` (defined by their derivative, improper
  by construction). Normalizing constants are kept throughout, so a
  proper penalty is exactly the negative log-density of its prior.
  `check_penalty()` validates every closed form against a route sharing
  no code with it.
