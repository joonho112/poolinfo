# Posterior density of a diagnostic

Computes the exact posterior density of a diagnostic or of the
between-site standard deviation on a grid of values.

## Usage

``` r
functional_density(
  x,
  which = c("I", "mean", "H", "J_eff", "sigma_tau"),
  n = 512
)
```

## Arguments

- x:

  An object of class `reliability_posterior`.

- which:

  The quantity: one of `"I"` (the default), `"mean"`, `"H"`, `"J_eff"`
  or `"sigma_tau"`.

- n:

  Number of points of the grid, an integer of at least 10.

## Value

A data frame with the columns `value` and `density`.

## Details

The density follows from that of `sigma_tau` by a change of variables:
for an increasing quantity \\f\\, the density at \\f(\sigma)\\ is
\\g(\sigma) / f'(\sigma)\\. For the dispersion \\H\\ the terms of the
two branches on either side of its largest value are added, with
\\dH/d\sigma = 2(\bar S - I)/\sigma\\. Under every prior of
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md)
the posterior density of `sigma_tau` has a finite positive limit at 0.
The densities of \\I\\ and the mean reliability grow without bound at
their values for `sigma_tau = 0`, like the inverse square root of the
distance from that value; they are still integrable. This behavior also
applies to \\H\\ and the concentration count when the standard errors
are unequal. With equal standard errors, their posteriors are point
masses at 0 and \\J\\, respectively.

The grid spans the central 99.9% of the posterior, so a plot of the
density does not include the extreme tails. With equal sampling
variances, \\J\_{\mathrm{eff}}\\ is a point mass at \\J\\, has no
density, and asking for one is an error. When the sampling variances
differ so little (by about 1e-6 relative or less) that the derivative of
\\J\_{\mathrm{eff}}\\ cannot be computed reliably, asking for its
density is an error too.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
head(functional_density(fit, "I", n = 50))
#>          value     density
#> 1 1.599263e-07 1563.220044
#> 2 1.913529e-02    4.543663
#> 3 3.827043e-02    3.226804
#> 4 5.740556e-02    2.643285
#> 5 7.654069e-02    2.294170
#> 6 9.567583e-02    2.054244
```
