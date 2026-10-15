# Site table of a profile

Returns one row for each site: its estimate and standard error, its
reliability, its term of the dispersion \\H\\, its proportion of the
total reliability, and, when the estimates are known, the posterior mean
and standard deviation of its effect.

## Usage

``` r
site_table(x, level = NULL, ...)
```

## Arguments

- x:

  An object of class `reliability_profile`, `reliability_posterior`,
  `reliability_draws` or `reliability_design`.

- level:

  Probability content of the posterior intervals, a number between 0
  and 1. `NULL`, the default, uses the level of the fit. Not used for a
  profile at a fixed between-site SD.

- ...:

  Must be empty.

## Value

A data frame with one row per site and the columns

- `site`, `estimate` (when known) and `se`.

- for a design, `n`, `kappa` and `n_balanced`: the planned size, the
  variance constant and the size under balance or equal spending.

- `reliability`: the site's reliability, its posterior mean for a
  posterior fit or draws, which also have `reliability_lower` and
  `reliability_upper`, the ends of its posterior interval.

- `kl_term`: the site's term of \\H\\; `share`: its proportion of the
  total reliability.

- when the estimates are known, `effect` and `effect_sd`: the posterior
  mean and standard deviation of the site's effect, and for a posterior
  fit or draws `effect_lower` and `effect_upper`, the ends of its
  posterior interval.

## Details

For a profile at a fixed `sigma_tau`, the effect's posterior mean and
standard deviation are those of the normal model at that `sigma_tau`,
with a flat prior on the mean of the effects: \\S_j \hat\tau_j + (1 -
S_j)\hat\mu\\ and \\\\\sigma\_\tau^2(1 - S_j) + (1 - S_j)^2 /
W\\^{1/2}\\, where \\\hat\mu\\ is the precision-weighted mean of the
estimates with weights \\w_j = 1/(\sigma\_\tau^2 + v_j)\\ and \\W =
\sum_j w_j\\. A profile inherited from a fitted mean model uses the
site's fitted value \\x_j'\hat\beta\\ in place of the common mean, and
\\x_j'(X'\Omega X)^{-1}x_j\\ in place of \\1/W\\, where \\\Omega\\ is
the diagonal matrix of the weights \\w_j\\.

The column `share` is the proportion of the total reliability at site
\\j\\, \\S_j / \sum_k S_k\\. At a fixed `sigma_tau` the inverse of the
sum of the squared proportions is the concentration count. For a
posterior fit the columns are posterior means, and that identity is only
approximate.

## See also

[`profile_summary()`](https://joonho112.github.io/poolinfo/reference/profile_summary.md)
for the diagnostics, and
[`plot_pooling()`](https://joonho112.github.io/poolinfo/reference/plot_pooling.md)
and
[`plot_profile()`](https://joonho112.github.io/poolinfo/reference/plot_profile.md)
for plots of the sites.

## Examples

``` r
prof <- reliability_profile(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2,
                            estimate = c(0.3, 0.1, -0.1, 0.5))
site_table(prof)
#>   site estimate   se reliability     kl_term     share     effect  effect_sd
#> 1    1      0.3 0.10   0.8000000 0.134361601 0.3559206 0.27629021 0.09333724
#> 2    2      0.1 0.15   0.6400000 0.009963631 0.2847365 0.12932238 0.12925311
#> 3    3     -0.1 0.20   0.5000000 0.010239076 0.2224504 0.04072553 0.15636178
#> 4    4      0.5 0.30   0.3076923 0.148205370 0.1368925 0.27946612 0.19032017
```
