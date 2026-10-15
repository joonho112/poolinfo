# Sites' terms of the dispersion

Plots each site's term \\\mathrm{KL}\\\mathrm{Ber}(I) \\
\mathrm{Ber}(S_j)\\\\, whose mean is the dispersion \\H\\, in order.

## Usage

``` r
plot_contributions(x, max_sites = 60, ...)
```

## Arguments

- x:

  An object of class `reliability_profile`, `reliability_posterior`,
  `reliability_draws` or `reliability_design`.

- max_sites:

  Largest number of sites to label, a positive integer.

- ...:

  Must be empty.

## Value

A ggplot object.

## Details

The terms are those of
[`kl_contributions()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md);
for a posterior fit they are posterior means.

## Examples

``` r
prof <- reliability_profile(se = eight_schools$se, sigma_tau = 10,
                            site = eight_schools$school)
plot_contributions(prof)
```
