# Plot of the posterior density of a diagnostic

Plots the exact posterior density of the informativeness, or of another
diagnostic, with its posterior mean, its value at the posterior mean of
`sigma_tau^2` and its central posterior interval.

## Usage

``` r
plot_functionals(
  x,
  which = c("I", "mean", "H", "J_eff", "sigma_tau"),
  level = NULL,
  plugin = TRUE,
  scale = 1,
  ...
)
```

## Arguments

- x:

  An object of class `reliability_posterior` or `reliability_draws`, or
  a `profile_comparison` of them.

- which:

  The diagnostic: one of `"I"` (the default), `"mean"`, `"H"`, `"J_eff"`
  or `"sigma_tau"`.

- level:

  Probability content of the shaded interval, a number between 0 and 1.
  `NULL`, the default, uses the level of the fit.

- plugin:

  If `TRUE` (the default), a dashed line marks the value at the
  posterior mean of `sigma_tau^2`, when that mean exists.

- scale:

  Factor that multiplies `sigma_tau` for `which = "sigma_tau"`, a
  positive number.

- ...:

  Must be empty.

## Value

A ggplot object.

## Details

For a posterior fit the density comes from
[`functional_density()`](https://joonho112.github.io/poolinfo/reference/functional_density.md).
The density of the informativeness grows without bound as \\I\\
approaches 0 under every prior offered by
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md),
whereas the posterior density of `sigma_tau` has a finite positive limit
at 0. For draws it is a kernel density estimate of the draws. For the
informativeness, a label gives the posterior probability that \\I\\
exceeds the threshold of the fit; the gray line is at one half, where
the second-order term of the mean reliability minus \\I\\ changes sign.
See
[`prob_above()`](https://joonho112.github.io/poolinfo/reference/prob_above.md).

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
plot_functionals(fit)

plot_functionals(fit, which = "H")
```
