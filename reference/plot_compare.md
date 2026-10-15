# Comparison plot

Plots the informativeness, the mean reliability and the dispersion of
several analyses as points with their posterior intervals, one row per
analysis and one panel per diagnostic.

## Usage

``` r
plot_compare(x, which = c("I", "mean", "H"), ...)
```

## Arguments

- x:

  An object of class `profile_comparison`.

- which:

  The diagnostics: any of `"I"`, `"mean"`, `"H"` and `"J_eff"`.

- ...:

  Must be empty.

## Value

A ggplot object.

## Details

The numbers are those of
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) of the
comparison. Analyses at a fixed `sigma_tau` have points without
intervals.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
cmp <- prior_sensitivity(fit, list(`half-normal` = sd_prior("half_normal", scale = 10)))
plot_compare(cmp)
```
