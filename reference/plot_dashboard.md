# Dashboard of a posterior fit

Combines four plots of a posterior fit in the layout of Figure 8 of Lee
(2026): (a) the pooling plot, (b) the reliability profile, (c) the
posterior density of the informativeness and (d) the informativeness and
logit variance plane.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a posterior
fit shows it.

## Usage

``` r
plot_dashboard(x, level = NULL, max_sites = 60, scale = 1, ...)
```

## Arguments

- x:

  An object of class `reliability_posterior` or `reliability_draws`.

- level:

  Probability content of the intervals, a number between 0 and 1.
  `NULL`, the default, uses the level of the fit.

- max_sites:

  Largest number of sites to label, a positive integer.

- scale:

  Factor that multiplies the estimates and the effects in panel (a), a
  positive number.

- ...:

  Must be empty.

## Value

A patchwork object.

## Details

Draws without estimates or without draws of the mean have no panel (a).
The captions of the panels are gathered below the plot.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
plot_dashboard(fit)
```
