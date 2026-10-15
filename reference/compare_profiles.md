# Diagnostics of several analyses side by side

Collects several posterior fits, sets of draws, profiles or designs,
such as the analyses of several outcomes, effect measures or subsets of
sites, and tabulates their diagnostics side by side.

## Usage

``` r
compare_profiles(..., .list = NULL)

# S3 method for class 'profile_comparison'
plot(x, type = NULL, ...)
```

## Arguments

- ...:

  For `compare_profiles()`, named objects of class
  `reliability_posterior`, `reliability_draws`, `reliability_profile` or
  `reliability_design`. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), further
  arguments of the plot function.

- .list:

  Optionally, a named list of such objects, added after those in `...`.

- x:

  A comparison made by `compare_profiles()`.

- type:

  The plot of a comparison: `"functionals"` (the posterior densities of
  the informativeness, for posterior fits and draws only), `"plane"`,
  `"compare"`, `"profile"`, `"gap"` or `"concentration"`. `NULL`, the
  default, gives `"functionals"` when every analysis is a posterior fit
  or draws, and `"compare"` otherwise.

## Value

An object of class `profile_comparison`, a list with the components

- `objects`: the named list of analyses.

- `table`: a data frame with one row per analysis and the columns
  `analysis`, `source`, `J`, `I`, `I_lower`, `I_upper`, `mean`,
  `mean_lower`, `mean_upper`, `H`, `H_lower`, `H_upper`, `J_eff`, `V`,
  `prob_I_above` and `sigma_tau` (posterior means, or values at a fixed
  `sigma_tau`, with the ends of the posterior intervals).

[`print()`](https://rdrr.io/r/base/print.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (the
table),
[`profile_summary()`](https://joonho112.github.io/poolinfo/reference/profile_summary.md),
[`report_table()`](https://joonho112.github.io/poolinfo/reference/report_table.md),
[`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) take it.

## Details

The table has one row per analysis with the informativeness, the mean
reliability, the dispersion, the concentration count and the logit
variance, and their posterior intervals where there are any. The
analyses can be of different kinds, which the column `source` records.
Each analysis is conditional on its own standard errors, so analyses of
the same sites with different effect measures (a risk difference and a
log odds ratio, for example) have different diagnostics.

## See also

[`prior_sensitivity()`](https://joonho112.github.io/poolinfo/reference/prior_sensitivity.md)
for the same analysis under several priors, and
[`vignette("multisite-trial")`](https://joonho112.github.io/poolinfo/articles/multisite-trial.md)
for a comparison of the analyses of one trial.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
prof <- reliability_profile(fit)
cmp <- compare_profiles(posterior = fit, `plug-in` = prof)
cmp
#> Comparison of 2 analyses
#> Intervals are central posterior intervals.
#> Each analysis is conditional on its own standard errors.
#> 
#>   analysis   J  informativeness I     mean reliability
#>   posterior  8  0.225 [0.000, 0.745]  0.230
#>   plug-in    8  0.338                 0.346
#> 
#>   analysis   dispersion H             concentration count J_eff
#>   posterior  0.0143 [0.0000, 0.0283]  7.2
#>   plug-in    0.0249                   7.4
```
