# Draws from the posterior

Draws the between-site standard deviation, the mean of the effects and
the diagnostics from their posterior, and optionally the site effects.

## Usage

``` r
profile_draws(x, ndraws = 4000, seed = NULL, sites = FALSE)
```

## Arguments

- x:

  An object of class `reliability_posterior`.

- ndraws:

  Number of draws, a positive integer.

- seed:

  Seed of the random number generator, an integer, or `NULL` (the
  default) for a seed from the clock. The seed used is stored in the
  attribute `seed`. The state of the session's generator is restored
  afterwards, so [`set.seed()`](https://rdrr.io/r/base/Random.html)
  before the call does not make the draws reproducible; give `seed`
  instead.

- sites:

  If `TRUE`, the draws include the site effects.

## Value

A data frame with one row per draw and the columns `draw`, `sigma_tau`,
`effect_mean` (or one column per coefficient of the mean model), `I`,
`mean`, `H`, `J_eff` and, with `sites = TRUE`, one column per site,
named `effect_` followed by the site's label (`effect_A` for site A).

## Details

`sigma_tau` is drawn by inverting a table of its distribution function,
then the mean of the effects from its normal distribution given
`sigma_tau`, and the site effects from their normal distributions given
both (Lee, 2026, Appendix D.3). The draws are independent. They are for
plots and for summaries that the package does not compute; the summaries
of the package itself come from the integrals and have no Monte Carlo
error.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
draws <- profile_draws(fit, ndraws = 1000, seed = 1)
mean(draws$I)
#> [1] 0.2258355
```
