# Three centroids of a reliability profile

Computes the harmonic centroid, the informativeness and the mean
reliability of a set of site reliabilities.

## Usage

``` r
reliability_centroids(S)
```

## Arguments

- S:

  Site reliabilities, a numeric vector of numbers strictly between 0 and
  1, as returned by
  [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md).

## Value

A named numeric vector with elements `harmonic`, `I` (the
informativeness) and `mean`.

## Details

The three numbers are the first-argument centroids of the reliabilities
under three Bregman divergences (Lee, 2026, Table 2): the harmonic mean
(Burg potential), the informativeness \\I\\ (Bernoulli Kullback-Leibler)
and the arithmetic mean \\\bar S\\ (squared error). The harmonic
centroid is never larger than the other two, and equals them only when
all the reliabilities are equal. There is no general order between \\I\\
and the mean reliability: either can be the larger, depending on the
level and the skewness of the profile, and the sign of their difference
is reported by
[`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md).

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`bregman_centroid()`](https://joonho112.github.io/poolinfo/reference/bregman_centroid.md)
for the centroids on either side of each divergence.

## Examples

``` r
S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
reliability_centroids(S)
#>  harmonic         I      mean 
#> 0.4961240 0.5714286 0.5619231 
```
