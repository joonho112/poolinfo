# Dispersion and centroid gap along the informativeness

Plots the dispersion \\H\\ and the mean reliability minus \\I\\ as the
informativeness moves with the centered logits of the profile held
fixed, as they do when only `sigma_tau` changes.

## Usage

``` r
plot_gap_path(x, leading = TRUE, ...)
```

## Arguments

- x:

  An object of class `reliability_profile`, `reliability_posterior`,
  `reliability_draws` or `reliability_design`, or a `profile_comparison`
  of them.

- leading:

  If `TRUE` (the default), dashed gray lines show the second-order terms
  of Proposition 2.

- ...:

  Must be empty.

## Value

A patchwork object with two panels.

## Details

Along the path \\dH/dm = \bar S - I\\, where \\m\\ is the logit of \\I\\
(Lee, 2026, Proposition 2). For unequal sampling variances, \\H\\ is
strictly concave in \\I\\ along the path and is largest at the single
point where the gap crosses 0. With equal sampling variances, \\H\\ and
the gap are 0 throughout the path and there is no turning point. The
point on each curve marks the profile's own informativeness (for a
posterior fit, its posterior median), and the cross marks the turning
point.

## Examples

``` r
prof <- reliability_profile(se = eight_schools$se, sigma_tau = 10)
plot_gap_path(prof)
```
