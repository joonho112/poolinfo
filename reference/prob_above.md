# Posterior probability that a diagnostic exceeds a value

Computes the posterior probability that the informativeness, the mean
reliability, the concentration count or the between-site standard
deviation is above a given value.

## Usage

``` r
prob_above(x, threshold = NULL, which = c("I", "mean", "J_eff", "sigma_tau"))
```

## Arguments

- x:

  An object of class `reliability_posterior` or `reliability_draws`.

- threshold:

  The value, a single number. For `which = "I"`, `NULL` (the default)
  uses the `threshold` of the fit; for the other quantities a value must
  be given.

- which:

  The quantity: one of `"I"` (the default), `"mean"`, `"J_eff"` or
  `"sigma_tau"`.

## Value

A number between 0 and 1. For draws, its attribute `mcse` is the Monte
Carlo standard error.

## Details

The four quantities increase with `sigma_tau`, so the probability is
\\1 - F(\sigma_c)\\, where \\F\\ is the posterior distribution function
of `sigma_tau` and \\\sigma_c\\ the value at which the quantity equals
the threshold. For the informativeness and a threshold \\c\\, \\\sigma_c
= \\\mathrm{GM}(v)\\ c/(1-c)\\^{1/2}\\. When all sampling variances are
equal, the concentration count is the constant \\J\\; its probability of
strictly exceeding \\c\\ is 1 for \\c \< J\\ and 0 for \\c \ge J\\.

Lee (2026) reports the probability that \\I\\ exceeds 1/2, the value at
which the between-site variance equals the geometric mean of the
sampling variances. Above it the leading term of the mean reliability
minus \\I\\ is negative (Proposition 2), although the remainder can
change the sign of the difference. The value 1/2 is a reference point,
not a standard that an analysis must meet.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`plot_functionals()`](https://joonho112.github.io/poolinfo/reference/plot_functionals.md),
whose label gives this probability.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
prob_above(fit)
#> [1] 0.1379141
prob_above(fit, threshold = 0.3)
#> [1] 0.3122938
prob_above(fit, threshold = 5, which = "sigma_tau")
#> [1] 0.5194775
```
