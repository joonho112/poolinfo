# Posterior summaries under several priors

Refits a posterior under each of several priors on the between-site
standard deviation and collects the fits in a comparison, to see how
much the diagnostics depend on the prior.

## Usage

``` r
prior_sensitivity(x, priors)
```

## Arguments

- x:

  An object of class `reliability_posterior`.

- priors:

  A named list of priors made by
  [`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md).

## Value

An object of class `profile_comparison` whose first analysis is `x`
itself, named after its prior.

## Details

With few sites the flat prior on the standard deviation can matter
(Gelman, 2006): its posterior has a heavy upper tail. A half-normal or
half-t prior whose scale is set on the scale of the effects is a common
alternative. The comparison shows how far the conclusions move; it does
not choose a prior. The scale of a proper prior has no default, because
it depends on the units of the estimates.

## References

Gelman, A. (2006). Prior distributions for variance parameters in
hierarchical models (comment on article by Browne and Draper). *Bayesian
Analysis, 1*(3), 515-534.
[doi:10.1214/06-BA117A](https://doi.org/10.1214/06-BA117A)

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
prior_sensitivity(fit, list(
  `half-normal, scale 10` = sd_prior("half_normal", scale = 10),
  `half-t, scale 10, 3 df` = sd_prior("half_t", scale = 10, df = 3)
))
#> Comparison of 3 analyses
#> Intervals are central posterior intervals.
#> Each analysis is conditional on its own standard errors.
#> 
#>   analysis                J  informativeness I     mean reliability
#>   flat                    8  0.225 [0.000, 0.745]  0.230
#>   half-normal, scale 10   8  0.157 [0.000, 0.570]  0.162
#>   half-t, scale 10, 3 df  8  0.156 [0.000, 0.587]  0.161
#> 
#>   analysis                dispersion H             concentration count J_eff
#>   flat                    0.0143 [0.0000, 0.0283]  7.2
#>   half-normal, scale 10   0.0118 [0.0000, 0.0282]  7.1
#>   half-t, scale 10, 3 df  0.0116 [0.0000, 0.0282]  7.1
```
