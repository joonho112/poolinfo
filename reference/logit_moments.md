# Logit variance and third absolute moment

Computes the variance \\V\\ and the third absolute moment \\\rho_3\\ of
the logits of the site reliabilities.

## Usage

``` r
logit_moments(S)
```

## Arguments

- S:

  Site reliabilities, a numeric vector of numbers strictly between 0 and
  1, as returned by
  [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md).

## Value

A named numeric vector with elements `V` and `rho3`.

## Details

With \\\Delta_j = \mathrm{logit}\\ S_j - \mathrm{logit}\\ I\\, \$\$V =
\frac{1}{J}\sum_j \Delta_j^2, \qquad \rho_3 = \frac{1}{J}\sum_j
\|\Delta_j\|^3.\$\$ Both divide by \\J\\. Because \\\mathrm{logit}\\ S_j
= \log \sigma\_\tau^2 - \log v_j\\, the deviations, and so \\V\\ and
\\\rho_3\\, depend on the sampling variances only: they do not change
with the between-site standard deviation.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
logit_moments(S)
#>         V      rho3 
#> 0.6448550 0.6748889 

# The same at another sigma_tau: the logit moments do not change.
logit_moments(site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 1))
#>         V      rho3 
#> 0.6448550 0.6748889 
```
