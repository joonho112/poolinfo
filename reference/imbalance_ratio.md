# Ratio of the geometric to the arithmetic mean site size

Gives the ratio of the geometric to the arithmetic mean of the site
sizes, the factor by which unequal sizes shrink the size that enters the
informativeness, as a function of their coefficient of variation.

## Usage

``` r
imbalance_ratio(cv, type = c("lognormal", "local"))
```

## Arguments

- cv:

  Coefficients of variation of the site sizes, numbers of at least 0.

- type:

  `"lognormal"` (the default) for the exact ratio in a population of
  sizes (a lognormal population), \\(1 + \mathrm{CV}^2)^{-1/2}\\;
  `"local"` for the approximation \\\exp(-\mathrm{CV}^2/2)\\ of a finite
  set of sizes close to their mean.

## Value

A numeric vector of ratios, one per value of `cv`.

## Details

The two are different statements (Lee, 2026, Appendix E.2). The
lognormal ratio is an identity for a population; it does not give the
geometric mean of a finite set of sizes from their coefficient of
variation. The local approximation applies to a finite set, with an
explicit remainder, when the sizes are close to their mean, and needs
the coefficient of variation with divisor \\J\\, not \\J - 1\\.
[`reliability_design()`](https://joonho112.github.io/poolinfo/reference/reliability_design.md)
reports the exact ratio of a planned allocation. At a coefficient of
variation of 0.75 the lognormal ratio is 0.80 and the local
approximation 0.755.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
imbalance_ratio(c(0.25, 0.5, 0.75))
#> [1] 0.9701425 0.8944272 0.8000000
imbalance_ratio(c(0.25, 0.5, 0.75), type = "local")
#> [1] 0.9692332 0.8824969 0.7548396
```
