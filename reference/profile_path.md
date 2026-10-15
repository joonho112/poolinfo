# Dispersion and centroid gap along a path of informativeness

Keeps the centered logits of a reliability profile fixed and moves its
informativeness \\I\\, as a change of the between-site standard
deviation does when the sampling variances stay fixed. Returns the
dispersion, the mean reliability and their second-order terms at each
point of the path.

## Usage

``` r
profile_path(x, I = NULL, n = 201)
```

## Arguments

- x:

  Site reliabilities, a numeric vector of numbers strictly between 0 and
  1, or an object of class `reliability_profile`,
  `reliability_posterior`, `reliability_draws` or `reliability_design`.

- I:

  Values of the informativeness at which to evaluate the path, a numeric
  vector of numbers strictly between 0 and 1. By default `n` equally
  spaced values.

- n:

  Number of points of the default path, an integer of at least 2.

## Value

A data frame with one row for each value of `I` and the columns

- `I` and `logit_I`.

- `sigma_tau`, for an object.

- `H`, `mean` and `gap`: the dispersion, the mean reliability and the
  mean reliability minus `I`.

- `H_leading`, `gap_leading`, `H_bound` and `gap_bound`: the
  second-order terms and the bounds on the remainders, as in
  [`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md).

## Details

Write \\m\\ for the logit of \\I\\ and \\\Delta_j\\ for the centered
logits of the profile. Along the path \\\Delta_j\\ stay fixed and \\m\\
moves, and Proposition 2 of Lee (2026) gives \\\partial H / \partial m =
\bar S - I\\. So the dispersion rises from 0 while the mean reliability
is above \\I\\, falls back to 0 after it, and is largest where the gap
changes sign.

For an object, the path also gives the between-site standard deviation
at each point, \\\sigma\_\tau = \\\mathrm{GM}(v)\\ I/(1-I)\\^{1/2}\\,
where \\\mathrm{GM}(v)\\ is the geometric mean of the sampling
variances. For reliabilities alone the scale is unknown, and that column
is absent.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
path <- profile_path(S, n = 9)
path[, c("I", "H", "gap")]
#>     I          H           gap
#> 1 0.1 0.03020018  2.291741e-02
#> 2 0.2 0.05165415  2.834731e-02
#> 3 0.3 0.06605218  2.359490e-02
#> 4 0.4 0.07435859  1.310160e-02
#> 5 0.5 0.07707534 -1.387779e-17
#> 6 0.6 0.07435859 -1.310160e-02
#> 7 0.7 0.06605218 -2.359490e-02
#> 8 0.8 0.05165415 -2.834731e-02
#> 9 0.9 0.03020018 -2.291741e-02
```
