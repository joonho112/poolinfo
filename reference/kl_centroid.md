# Informativeness of a reliability profile

Computes the informativeness \\I\\, the centroid of the site
reliabilities under the Bernoulli Kullback-Leibler divergence: the
logistic function of the mean logit of the reliabilities.

## Usage

``` r
kl_centroid(S)
```

## Arguments

- S:

  Site reliabilities, a numeric vector of numbers strictly between 0 and
  1, as returned by
  [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md).

## Value

The informativeness, a number between 0 and 1.

## Details

The informativeness is \$\$I = \mathrm{logistic}\Big(\frac{1}{J}\sum_j
\mathrm{logit}\\ S_j\Big),\$\$ the value \\p\\ that minimizes the
average divergence \\J^{-1}\sum_j \mathrm{KL}\\\mathrm{Ber}(p) \\
\mathrm{Ber}(S_j)\\\\ (Lee, 2026, Theorem 1). When the reliabilities
come from sampling variances \\v_j\\ and a between-site standard
deviation \\\sigma\_\tau\\, it equals \\\sigma\_\tau^2 /
\\\sigma\_\tau^2 + \mathrm{GM}(v)\\\\, where \\\mathrm{GM}(v)\\ is the
geometric mean of the variances.

A reliability of exactly 0 or 1 has no finite logit, and such a profile
has no centroid; the function then fails with an error.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`reliability_centroids()`](https://joonho112.github.io/poolinfo/reference/reliability_centroids.md)
for the harmonic and arithmetic means beside \\I\\, and
[`kl_dispersion()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md)
for the spread around it.

## Examples

``` r
S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
kl_centroid(S)
#> [1] 0.5714286
```
