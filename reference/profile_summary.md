# Summaries of a profile as a data frame

Returns the diagnostics of a reliability profile, or their posterior
means and intervals, as a data frame with one row for each quantity.

## Usage

``` r
profile_summary(x, level = NULL, ...)
```

## Arguments

- x:

  An object of class `reliability_profile`, `reliability_posterior`,
  `reliability_draws`, `reliability_design` or `profile_comparison`.

- level:

  Probability content of the posterior intervals, a number between 0
  and 1. `NULL`, the default, uses the level of the fit. Not used for a
  profile at a fixed between-site SD.

- ...:

  Must be empty.

## Value

A data frame with one row per quantity and the columns

- `quantity`: the code of the quantity.

- `mean`, `sd`, `lower`, `upper`: the posterior mean and standard
  deviation and the ends of the central posterior interval.

- for draws, `mcse` and `ess`: the Monte Carlo standard error of the
  mean and the effective sample size of the draws.

- for a design, `allocation`: the planned allocation, or the balanced
  one (equal spending with costs).

- for a comparison, `analysis`, the name of the analysis, and `mcse`.

## Details

The quantities are coded as in the rest of the package: `I`
(informativeness), `mean` (mean reliability), `harmonic` (harmonic
centroid), `H` (dispersion), `V` (logit variance), `J_eff`
(concentration count) and `sigma_tau` (between-site standard deviation).
A posterior fit or draws add `sigma_tau2` and either `effect_mean`, the
mean of the effects, or `coef_<name>` for each coefficient of a mean
model. For a profile at a fixed `sigma_tau`, the column `mean` contains
the value of each quantity and the other columns are `NA`.

## See also

[`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md)
for the sites, and
[`report_table()`](https://joonho112.github.io/poolinfo/reference/report_table.md)
for a table for a report.

## Examples

``` r
prof <- reliability_profile(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
profile_summary(prof)
#>    quantity       mean sd lower upper
#> 1         I 0.57142857 NA    NA    NA
#> 2      mean 0.56192308 NA    NA    NA
#> 3  harmonic 0.49612403 NA    NA    NA
#> 4         H 0.07569242 NA    NA    NA
#> 5         V 0.64485497 NA    NA    NA
#> 6     J_eff 3.62347623 NA    NA    NA
#> 7 sigma_tau 0.20000000 NA    NA    NA
```
