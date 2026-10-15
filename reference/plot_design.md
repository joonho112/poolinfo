# Design plots

Plots what a planned allocation implies: the cost of unequal site sizes,
the trade between the number of sites and the precision of the mean at a
fixed total, or the site size needed for a target informativeness.

## Usage

``` r
plot_design(x, type = c("imbalance", "sites", "target"), ...)
```

## Arguments

- x:

  An object of class `reliability_design`.

- type:

  The plot: `"imbalance"` (the default), `"sites"` or `"target"`. See
  Details.

- ...:

  Must be empty.

## Value

A ggplot object.

## Details

- `"imbalance"`: the ratio of the geometric to the arithmetic mean of
  the sizes against their coefficient of variation, for a lognormal
  population of sizes (exact) and as the local approximation
  \\\exp(-\mathrm{CV}^2/2)\\, with the planned allocation marked (Lee,
  2026, Appendix E.2).

- `"sites"`: at the planned total \\N\\, the informativeness of a
  balanced design and the standard error of the mean effect,
  \\(\sigma\_\tau^2/J + \kappa/N)^{1/2}\\, against the number of sites
  \\J\\ (Appendix E.3).

- `"target"`: the informativeness of a balanced design against the size
  per site, one curve for each `sigma_tau` of the design.

The curves use the geometric mean of the design's \\\kappa_j\\.

## Examples

``` r
n <- round(stats::qlnorm(stats::ppoints(40), 0, 0.47) / 40 * 2000)
d <- reliability_design(n, sigma_tau = c(0.15, 0.2, 0.25), kappa = 4)
plot_design(d)

plot_design(d, type = "sites")

plot_design(d, type = "target")
```
