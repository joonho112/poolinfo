# Confidence interval for the between-site variance mapped through the diagnostics

Maps a confidence interval for the between-site variance, such as the
Q-profile interval of a metafor fit, through the informativeness, the
mean reliability, the concentration count and the dispersion. In
metafor's notation the between-site variance is \\\tau^2\\.

## Usage

``` r
map_tau2_interval(
  x = NULL,
  se = NULL,
  se2 = NULL,
  interval = NULL,
  level = 0.95
)
```

## Arguments

- x:

  An optional `rma.uni` fit of the metafor package, from which the
  sampling variances, the estimate of \\\tau^2\\ and its Q-profile
  confidence interval are computed.

- se, se2:

  Standard errors or sampling variances of the site estimates, when `x`
  is not given. Give one of them.

- interval:

  The confidence interval for the between-site variance, two numbers of
  at least 0, optionally with the point estimate, as in
  `c(estimate = 0.02, lower = 0.005, upper = 0.08)`. It is needed
  without `x`; with `x` it replaces the computed Q-profile interval.

- level:

  The level of the Q-profile interval computed for `x`, a number between
  0 and 1, when `interval` is not given.

## Value

A data frame with one row for each of `"I"`, `"mean"`, `"J_eff"` and
`"H"`, and the columns `functional`, `estimate` (at the point estimate
of \\\tau^2\\, when it is known), `lower` and `upper`. Its attribute
`interval` is the interval for \\\tau^2\\ that was mapped, and its
attribute `method` describes where that interval came from.

## Details

Confidence-interval mapping is an experimental extension beyond Lee
(2026); its interface may change in a later 0.x release.

The informativeness and the mean reliability increase with \\\tau^2\\,
so the interval between their values at the ends of a confidence
interval \\C\\ has the coverage of \\C\\. So does the interval of the
concentration count when the sampling variances differ; when they are
all equal the count is \\J\\ whatever \\\tau^2\\ is, and its interval is
the single value \\J\\. The dispersion \\H\\ first rises and then falls;
its interval is the image \\\[\min_C H, \max_C H\]\\, which covers the
true value at least as often as \\C\\ covers \\\tau^2\\. These are
frequentist intervals, conditional on the sampling variances like the
rest of the package; they are not posterior intervals.

For a metafor fit the interval is the Q-profile interval for the model
of the fit, moderators included. poolinfo computes it for every
estimator, also for `GENQ` and `GENQM` fits, whose own
[`confint()`](https://rdrr.io/r/stats/confint.html) uses another method,
and the result does not depend on the units of the effects. A fit
without a between-site variance (`method = "FE"`, `"EE"` or `"CE"`) is
an error.

## See also

[metafor-methods](https://joonho112.github.io/poolinfo/reference/metafor-methods.md)
for profiles and posteriors of metafor fits, and
[`vignette("meta-analysis")`](https://joonho112.github.io/poolinfo/articles/meta-analysis.md)
for an example.

## Examples

``` r
map_tau2_interval(se = eight_schools$se, interval = c(estimate = 20, lower = 0, upper = 500))
#>   functional   estimate    lower      upper
#> 1          I 0.11956387 0.000000 0.77246965
#> 2       mean 0.12840963 0.000000 0.76192310
#> 3      J_eff 7.02118867 6.817468 7.89281710
#> 4          H 0.01160368 0.000000 0.02836966
```
