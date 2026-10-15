# Sizes that reach a target informativeness

Computes, for a balanced design, the site size needed for a target
informativeness with a given number of sites, or the largest number of
sites that reaches the target with a given total.

## Usage

``` r
solve_design(target_I, sigma_tau, kappa, J = NULL, N = NULL)
```

## Arguments

- target_I:

  Target informativeness, a number strictly between 0 and 1.

- sigma_tau:

  Planned between-site standard deviation, a positive number, or a
  numeric vector of them for a range of scenarios.

- kappa:

  Constant \\\kappa\\ of the sampling variances \\v_j = \kappa / n_j\\,
  a positive number (see
  [`neyman_variance()`](https://joonho112.github.io/poolinfo/reference/neyman_variance.md)).

- J, N:

  The number of sites or the total number of individuals: give one of
  them.

## Value

A data frame with one row per value of `sigma_tau` and the columns

- `target_I`, `sigma_tau`, `kappa`: the inputs.

- `J`, `n_per_site`, `N`: the design, with `n_per_site` not rounded.

- `I`, with `N` given: the informativeness reached with the number of
  sites rounded down. `J`, `n_per_site` and `I` are `NA` when the total
  cannot reach the target even in a single site.

## Details

Under balance \\I = \sigma\_\tau^2 / (\sigma\_\tau^2 + \kappa / \bar
n)\\, so the target needs \\\bar n = (\kappa / \sigma\_\tau^2)\\ I/(1 -
I)\\ individuals per site, and with a total \\N\\ at most \\J = N
\sigma\_\tau^2 (1 - I) / (\kappa I)\\ sites, rounded down. Under balance
the informativeness depends on the size per site, not on the number of
sites: with `J` given, the target fixes the size per site, and with `N`
given, it fixes how many sites the total can be divided into. Because
`sigma_tau` is unknown at the planning stage, a vector of plausible
values gives one row per scenario.

At a fixed total, more sites lower the informativeness but make the mean
effect more precise: its variance under balance is \\\sigma\_\tau^2/J +
\kappa/N\\ (Appendix E.3). The two are different design questions.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`reliability_design()`](https://joonho112.github.io/poolinfo/reference/reliability_design.md)
and
[`plot_design()`](https://joonho112.github.io/poolinfo/reference/plot_design.md),
and
[`vignette("planning")`](https://joonho112.github.io/poolinfo/articles/planning.md).

## Examples

``` r
solve_design(target_I = 0.5, sigma_tau = c(0.1, 0.2, 0.3), kappa = 4, J = 40)
#>   target_I sigma_tau kappa  J n_per_site         N
#> 1      0.5       0.1     4 40  400.00000 16000.000
#> 2      0.5       0.2     4 40  100.00000  4000.000
#> 3      0.5       0.3     4 40   44.44444  1777.778
solve_design(target_I = 0.5, sigma_tau = 0.2, kappa = 4, N = 4000)
#>   target_I sigma_tau kappa  J n_per_site    N   I
#> 1      0.5       0.2     4 40        100 4000 0.5
```
