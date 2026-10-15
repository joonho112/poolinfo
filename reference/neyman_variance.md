# Sampling variances of a planned trial

Computes the planning approximation of the sampling variance of each
site's impact estimate in an individually randomized multisite trial.

## Usage

``` r
neyman_variance(n, p = 0.5, R2 = 0, var_outcome = 1)
```

## Arguments

- n:

  Planned numbers of individuals in the sites, a numeric vector of
  positive numbers.

- p:

  Proportion of each site's individuals assigned to treatment, numbers
  strictly between 0 and 1, one for all sites or one per site.

- R2:

  Proportion of the outcome variance explained by covariates, numbers in
  \\\[0, 1)\\, one for all sites or one per site.

- var_outcome:

  Variance of the outcome within each arm of a site (\\s_j^2\\ in
  equation 2 of Lee, 2026), positive numbers, one for all sites or one
  per site.

## Value

A numeric vector of sampling variances, one per site.

## Details

The variance of a difference in means with covariate adjustment is
approximately \$\$v_j = \Big(\frac{1}{p_j} + \frac{1}{1 -
p_j}\Big)\frac{s_j^2 (1 - R^2)}{n_j} = \frac{\kappa_j}{n_j}\$\$ (Lee,
2026, equation 2). The argument names are those of
[`multisiteDGP::compute_kappa()`](https://joonho112.github.io/multisiteDGP/reference/compute_kappa.html).
Designs with blocks, clusters or other estimators need their own
variance formula;
[`reliability_design()`](https://joonho112.github.io/poolinfo/reference/reliability_design.md)
also takes the constants \\\kappa_j\\ directly.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
neyman_variance(c(20, 50, 100))
#> [1] 0.20 0.08 0.04
neyman_variance(c(20, 50, 100), p = 0.6, R2 = 0.4)
#> [1] 0.125 0.050 0.025
```
