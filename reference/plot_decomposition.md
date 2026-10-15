# Decomposition of the average divergence

Plots the average divergence of a constant reliability \\p\\ from the
site reliabilities, which is smallest at \\p = I\\ with the value \\H\\,
and shows at one \\p\\ how it splits into \\\mathrm{KL}\\\mathrm{Ber}(p)
\\ \mathrm{Ber}(I)\\ + H\\.

## Usage

``` r
plot_decomposition(x, at = NULL, ...)
```

## Arguments

- x:

  An object of class `reliability_profile`, `reliability_posterior`,
  `reliability_draws` or `reliability_design`, or site reliabilities as
  a numeric vector.

- at:

  The constant reliability at which the split is shown, a number
  strictly between 0 and 1. By default a value well away from \\I\\.

- ...:

  Must be empty.

## Value

A ggplot object.

## Details

For every \\p\\, \\J^{-1}\sum_j
\mathrm{KL}\\\mathrm{Ber}(p)\\\mathrm{Ber}(S_j)\\ =
\mathrm{KL}\\\mathrm{Ber}(p)\\\mathrm{Ber}(I)\\ + H\\ (Lee, 2026,
Theorem 2). The mean reliability is marked too: it minimizes the average
divergence with \\p\\ in the second argument, \\J^{-1}\sum_j
\mathrm{KL}\\\mathrm{Ber}(S_j)\\\mathrm{Ber}(p)\\\\, not the average
plotted here. For a posterior fit the profile is taken at the posterior
median of `sigma_tau`.

## Examples

``` r
plot_decomposition(c(0.05, 0.10, 0.15, 0.25, 0.92))
```
