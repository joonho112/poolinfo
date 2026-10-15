# Profiles and posteriors of simulated multisite data

The data simulated by
[`multisiteDGP::sim_multisite()`](https://joonho112.github.io/multisiteDGP/reference/sim_multisite.html)
and
[`multisiteDGP::sim_meta()`](https://joonho112.github.io/multisiteDGP/reference/sim_meta.html)
carry the true between-site standard deviation of their design.
[`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
gives the profile at that true value, which an analysis of real data
could not know.
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
fits the posterior from the estimates and the standard errors alone, as
an analysis of real data would.

## Usage

``` r
# S3 method for class 'multisitedgp_data'
reliability_profile(x, ...)
```

## Arguments

- x:

  A `multisitedgp_data` object.

- ...:

  Must be empty.

## Value

A `reliability_profile` object whose source is the design of the
simulation.

## Details

The method reads the columns `site_index`, `tau_j_hat` and `se2_j` and
the element `sigma_tau` of the attribute `design`, which are documented
in
[`?multisiteDGP::multisitedgp_data`](https://joonho112.github.io/multisiteDGP/reference/multisitedgp_data.html).
The true effects `tau_j` are not used. Selecting rows retains the
simulation design. If the design attribute has been removed or has no
`sigma_tau`, supply the between-site standard deviation directly with
`sigma_tau`, together with `se2 = x$se2_j`.

## Examples

``` r
design <- multisiteDGP::multisitedgp_design(J = 30, sigma_tau = 0.2, seed = 1)
d <- multisiteDGP::sim_multisite(design)
reliability_profile(d)
#> Reliability profile of 30 sites at sigma_tau = 0.2 (multisiteDGP design value)
#> 
#>   informativeness I            0.301
#>   mean reliability             0.309
#>   dispersion H                 0.0189
#>   logit variance V             0.186
#>   concentration count J_eff    27.8 of 30
#> 
#> Site reliabilities range from 0.14 to 0.47.
reliability_posterior(d)
#> Reliability profile of 30 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#>                              mean     95% interval
#>   informativeness I          0.298    0.005    0.614
#>   mean reliability           0.303    0.006    0.610
#>   dispersion H               0.0164   0.0005   0.0228
#>   concentration count J_eff  27.7     26.1     29.2     of 30
#>   logit variance V           0.186    (fixed by the standard errors)
#>   sigma_tau                  0.200    0.0220   0.384
#>   effect mean                0.0339   -0.0979  0.171
#> 
#>   Pr(I > 0.5)                0.125
#>   I at E(sigma_tau^2)        0.343
```
