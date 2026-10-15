# Benchmark risk of a reliability profile

Computes \\\sigma\_\tau^2 (1 - \bar S)\\, the average squared error of
the best predictors that are affine in each site's own estimate when the
mean and the variances are known.

## Usage

``` r
benchmark_risk(..., se = NULL, se2 = NULL, sigma_tau)
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

The benchmark risk, a number of at least 0, in the squared units of the
estimates.

## Details

For a known common mean and known \\\sigma\_\tau^2\\ and \\v_j\\, the
predictor \\\mu + S_j(\hat\tau_j - \mu)\\ has the smallest mean squared
error among predictors affine in the site's own estimate, and its error
averaged over sites is \\\sigma\_\tau^2(1 - \bar S)\\ (Lee, 2026,
Proposition 2(iii)). Under the normal model with independent sites it is
the Bayes risk of the posterior means with known hyperparameters. It is
a benchmark at known inputs, not the error of a fitted model, whose
hyperparameters are estimated. It is 0 at `sigma_tau = 0` and approaches
the mean of the \\v_j\\ as `sigma_tau` grows.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
benchmark_risk(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#> [1] 0.01752308
```
