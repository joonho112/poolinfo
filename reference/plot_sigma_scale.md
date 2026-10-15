# Posterior of sigma_tau on the scale of the informativeness

Plots the posterior density of the between-site standard deviation with
an axis of the informativeness on top, and shades the region where the
informativeness exceeds the threshold of the fit.

## Usage

``` r
plot_sigma_scale(x, prior = TRUE, scale = 1, ...)
```

## Arguments

- x:

  An object of class `reliability_posterior` or `reliability_draws`.

- prior:

  If `TRUE` (the default), a proper prior is shown too, scaled to the
  height of the posterior.

- scale:

  Factor that multiplies `sigma_tau`, a positive number.

- ...:

  Must be empty.

## Value

A ggplot object.

## Details

The informativeness \\I = \sigma\_\tau^2/\\\sigma\_\tau^2 +
\mathrm{GM}(v)\\\\ increases with `sigma_tau`, so its posterior is the
image of that of `sigma_tau`; \\I \> c\\ exactly when \\\sigma\_\tau \>
\\\mathrm{GM}(v)\\ c/(1-c)\\^{1/2}\\, which is \\\mathrm{GM}(v)^{1/2}\\
for \\c = 1/2\\. For draws the density is a kernel estimate.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
plot_sigma_scale(fit)
```
