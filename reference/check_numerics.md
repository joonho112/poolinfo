# Numerical checks of a posterior fit

Reports the checks of the integration over the between-site standard
deviation: the agreement of the rule with two other rules, the bounds on
the omitted upper tail, and the location of the peak of the dispersion
found at two resolutions.

## Usage

``` r
check_numerics(x)
```

## Arguments

- x:

  An object of class `reliability_posterior`.

## Value

A data frame with the columns `check`, `value`, `limit` and `status`
(`"ok"` or `"outside"`).

## Details

For the node rule, the normalizing constant and the posterior means of
\\I\\ and \\H\\ are compared with those of a rule whose panels are all
halved, and with those of a rule whose lowest panel ends at half the
value. The differences are relative, and their limit is \\10^{-10}\\, or
\\64J\\ times the machine epsilon when that is larger, because the log
density is a sum of \\J\\ terms. For the adaptive engine, the check is
the error estimate of
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) for the
normalizing constant. The bounds on the omitted upper tail are those of
Appendix D.3 of Lee (2026), extended to the mean model and to the priors
of
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md);
the bound for \\E(\sigma\_\tau^q)\\ is shown divided by
\\\sqrt{\mathrm{GM}(v)}^q\\, so that it does not depend on the units.
The bounds concern the numerical approximation, not the calibration of
the posterior. The last row compares the peak of \\H\\, where the mean
reliability equals \\I\\, as located by scans of the logit of \\I\\ at
two resolutions.
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
has the details.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
check_numerics(fit)
#>                                                             check        value
#> 1     normalizing constant, rule against halved panels (relative) 0.000000e+00
#> 2            E(I) and E(H), rule against halved panels (relative) 1.214135e-16
#> 3           normalizing constant, lower end moved down (relative) 0.000000e+00
#> 4                  E(I) and E(H), lower end moved down (relative) 0.000000e+00
#> 5                               mass beyond the upper end (bound) 1.815325e-33
#> 6 E(sigma_tau^1) beyond the upper end, over sqrt(GM(v))^1 (bound) 7.171246e-28
#> 7 E(sigma_tau^2) beyond the upper end, over sqrt(GM(v))^2 (bound) 2.950961e-22
#> 8 E(sigma_tau^4) beyond the upper end, over sqrt(GM(v))^4 (bound) 6.396041e-11
#> 9            turning points of H at two resolutions (logit scale) 7.355228e-16
#>   limit status
#> 1 1e-10     ok
#> 2 1e-10     ok
#> 3 1e-10     ok
#> 4 1e-10     ok
#> 5 1e-12     ok
#> 6 1e-10     ok
#> 7 1e-10     ok
#> 8 1e-10     ok
#> 9 1e-07     ok
```
