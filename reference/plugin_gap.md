# Posterior means against plug-in values

Compares the posterior mean of each diagnostic with its value at the
posterior mean of \\\sigma\_\tau^2\\, the value that a plug-in of that
estimate would give.

## Usage

``` r
plugin_gap(x)
```

## Arguments

- x:

  An object of class `reliability_posterior` or `reliability_draws`.

## Value

A data frame with the columns `functional` (`"I"`, `"mean"`, `"H"` and
`"J_eff"`), `posterior_mean`, `at_mean_sigma2` and `difference` (the
second minus the first), and the attribute `sigma2`, the posterior mean
of \\\sigma\_\tau^2\\.

## Details

The informativeness and every site reliability are concave functions of
\\\sigma\_\tau^2\\. By Jensen's inequality, the posterior means of the
informativeness and of the mean reliability are therefore at most their
values at the posterior mean of \\\sigma\_\tau^2\\, and strictly less
unless the posterior of \\\sigma\_\tau^2\\ is a point mass (Lee, 2026,
Proposition 3). The dispersion and the concentration count have no such
general order. A plug-in of another estimate, such as the REML estimate,
can fall on either side of the posterior mean.

The comparison needs the posterior mean of \\\sigma\_\tau^2\\, which
under the flat prior exists only with more than 4 sites (more than \\p +
3\\ with a mean model of \\p\\ columns). Without it the values are `NA`,
with a message that gives the number of sites needed.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`vignette("posterior-summaries")`](https://joonho112.github.io/poolinfo/articles/posterior-summaries.md)
for the comparison in an example.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
plugin_gap(fit)
#>   functional posterior_mean at_mean_sigma2 difference
#> 1          I     0.22529451     0.33790826 0.11261375
#> 2       mean     0.22972844     0.34644286 0.11671441
#> 3          H     0.01428773     0.02491994 0.01063221
#> 4      J_eff     7.17004472     7.36307090 0.19302618
```
