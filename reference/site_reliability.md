# Site reliabilities

Computes the reliability of each site: the proportion of the variance of
its estimate that comes from the variation of the true effects across
sites. Under the normal hierarchical model with a known mean and a known
between-site standard deviation, it is also the weight of the site's own
estimate in its posterior mean.

## Usage

``` r
site_reliability(..., se = NULL, se2 = NULL, sigma_tau)
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

- sigma_tau:

  Standard deviation of the true effects across sites, a number of at
  least 0.

## Value

A numeric vector of reliabilities, with the names of `se` or `se2`.

## Details

The reliability of site \\j\\ is \$\$S_j =
\frac{\sigma\_\tau^2}{\sigma\_\tau^2 + v_j},\$\$ where \\v_j\\ is its
sampling variance. It equals
\\\mathrm{Var}(\tau_j)/\mathrm{Var}(\hat\tau_j)\\ whenever the estimate
is unbiased given the effect, without normality or independence across
sites (Lee, 2026, Proposition 1).

With `sigma_tau = 0` every reliability is 0. The functions of the
reliabilities, such as
[`kl_centroid()`](https://joonho112.github.io/poolinfo/reference/kl_centroid.md),
need them above 0, and
[`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
gives the limits of the diagnostics as `sigma_tau` goes to 0.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#> [1] 0.8000000 0.6400000 0.5000000 0.3076923
```
