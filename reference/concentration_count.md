# Concentration count of a reliability profile

Computes the concentration count \\J\_{\mathrm{eff}}\\, the inverse of
the sum of the squared proportions \\S_j / \sum_k S_k\\. It describes
how the total reliability is shared among the sites.

## Usage

``` r
concentration_count(S)
```

## Arguments

- S:

  Site reliabilities, a numeric vector of numbers strictly between 0 and
  1, as returned by
  [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md).

## Value

The concentration count, a number between 1 and the number of sites.

## Details

With \\q_j = S_j / \sum_k S_k\\, \$\$J\_{\mathrm{eff}} = \frac{1}{\sum_j
q_j^2} = \frac{(\sum_j S_j)^2}{\sum_j S_j^2}\$\$ (Lee, 2026, Section 5).
It lies between 1 and \\J\\, the number of sites, and equals \\J\\ only
when all the reliabilities are equal; it approaches 1 when one site has
nearly all of the reliability. It does not change when all the
reliabilities are multiplied by a common factor, so profiles that are
uniformly low and uniformly high both give \\J\\.

Its formula is Kish's, with the reliabilities as weights, but the count
is not a sample size or a number of degrees of freedom, and it does not
measure how well the sites can be ranked. A small value means that a few
sites carry much of the total reliability, and
[`plot_concentration()`](https://joonho112.github.io/poolinfo/reference/plot_concentration.md)
shows which.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
concentration_count(S)
#> [1] 3.623476

# One precise site and three noisy ones
concentration_count(c(0.95, 0.05, 0.05, 0.05))
#> [1] 1.32967
```
