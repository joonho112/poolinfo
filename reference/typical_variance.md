# Typical sampling variance of Higgins and Thompson

Computes the "typical" within-study variance of Higgins and Thompson
(2002), the variance that the \\I^2\\ statistic of meta-analysis
compares with the between-study variance.

## Usage

``` r
typical_variance(..., se = NULL, se2 = NULL)
```

## Arguments

- ...:

  Must be empty. The standard errors and the sampling variances are
  named arguments, because an unnamed vector could hold either.

- se:

  Standard errors of the site estimates, a numeric vector of positive
  numbers. Give `se` or `se2`, not both.

- se2:

  Sampling variances of the site estimates, a numeric vector of positive
  numbers.

## Value

The typical variance, a positive number.

## Details

With weights \\w_j = 1 / v_j\\, \\W = \sum_j w_j\\ and \\W_2 = \sum_j
w_j^2\\, \$\$v\_{\mathrm{typ}} = \frac{(J-1)W}{W^2 - W_2}.\$\$ For a
random-effects model without moderators, metafor's `rma()` reports
\\I^2\\ as the percentage \\100\\\hat\sigma\_\tau^2 /
(\hat\sigma\_\tau^2 + v\_{\mathrm{typ}})\\ at its estimate
\\\hat\sigma\_\tau^2\\ (\\\tau^2\\ in metafor's notation). The typical
variance is at least the harmonic mean of the \\v_j\\ (Böhning et al.,
2017), and in general it differs from their geometric mean, which
defines the informativeness. At least 2 sites are needed.

## References

Böhning, D., Lerdsuwansri, R., & Holling, H. (2017). Some general points
on the \\I^2\\-measure of heterogeneity in meta-analysis. *Metrika, 80*,
685-695.
[doi:10.1007/s00184-017-0622-3](https://doi.org/10.1007/s00184-017-0622-3)

Higgins, J. P. T., & Thompson, S. G. (2002). Quantifying heterogeneity
in a meta-analysis. *Statistics in Medicine, 21*, 1539-1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

## See also

[`kl_centroid()`](https://joonho112.github.io/poolinfo/reference/kl_centroid.md)
for the informativeness, which uses the geometric mean of the sampling
variances instead.

## Examples

``` r
se <- c(0.10, 0.15, 0.20, 0.30)
typical_variance(se = se)
#> [1] 0.02725155
exp(mean(log(se^2)))  # the geometric mean, for comparison
#> [1] 0.03
```
