# Concentration of the reliability

Plots the cumulative proportion of the total reliability against the
fraction of sites, the most reliable sites first, with the diagonal of
equal proportions.

## Usage

``` r
plot_concentration(x, label_top = 3, ...)
```

## Arguments

- x:

  An object of class `reliability_profile`, `reliability_posterior`,
  `reliability_draws` or `reliability_design`, or a `profile_comparison`
  of them.

- label_top:

  Number of the most reliable sites to name, an integer of at least 0.

- ...:

  Must be empty.

## Value

A ggplot object.

## Details

The proportions are \\q_j = S_j / \sum_k S_k\\, and for a fixed profile
the inverse of the sum of their squares is the concentration count
\\J\_{\mathrm{eff}}\\. For posterior fits and draws, each step adds the
posterior mean of a site's proportion, and the label gives the posterior
mean of \\J\_{\mathrm{eff}}\\, which in general differs from the inverse
of the sum of the squared posterior mean proportions. A curve on the
diagonal means that every site has the same proportion.

## Examples

``` r
prof <- reliability_profile(se = eight_schools$se, sigma_tau = 10,
                            site = eight_schools$school)
plot_concentration(prof)
```
