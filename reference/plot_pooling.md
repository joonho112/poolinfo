# Pooling plot

Plots each site's estimate with its interval beside the posterior mean
of its effect with its posterior interval, so that the partial pooling
of every site can be seen.

## Usage

``` r
plot_pooling(
  x,
  order = c("estimate", "reliability", "input"),
  level = NULL,
  color_by_reliability = FALSE,
  max_sites = 60,
  scale = 1,
  ...
)
```

## Arguments

- x:

  An object of class `reliability_posterior`, `reliability_draws` (with
  estimates and draws of the mean) or `reliability_profile` (with
  estimates).

- order:

  Order of the sites from bottom to top: `"estimate"` (the default),
  `"reliability"`, or `"input"` for the order of the input.

- level:

  Probability content of the intervals, a number between 0 and 1.
  `NULL`, the default, uses the level of the fit (0.95 for a profile).

- color_by_reliability:

  If `TRUE`, the posterior means are colored by the site reliability
  instead of a single color.

- max_sites:

  Largest number of sites to label, a positive integer. With more sites
  the axis has no labels.

- scale:

  Factor that multiplies the estimates and the effects, a positive
  number; `100` turns proportions into percentage points.

- ...:

  Must be empty.

## Value

A ggplot object.

## Details

The estimate's interval is \\\hat\tau_j \pm z\\ \mathrm{se}\_j\\, with
\\z\\ the standard normal quantile at \\(1 + \mathrm{level})/2\\. The
posterior intervals of a posterior fit are quantiles of the normal
mixture over `sigma_tau` of the effect's conditional posterior, so they
are not symmetric in general. For a profile at a fixed `sigma_tau`, the
posterior is normal, with the mean and standard deviation of
[`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md).
The vertical line is the posterior mean of the effect mean (for a
profile, the precision-weighted mean at `sigma_tau`). The horizontal
axis names the scale of the effects when the fit records it. With a mean
model there is no single mean: each posterior mean is pulled toward the
site's own predicted value, and there is no vertical line.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
plot_pooling(fit)
```
