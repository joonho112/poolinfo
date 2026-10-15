# Profiles and posteriors of metafor fits

[`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
of a random-effects fit of the metafor package gives the profile at its
estimate of \\\tau^2\\, a plug-in profile;
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
of the same fit integrates over the between-study standard deviation
instead, with the same estimates, sampling variances and moderators.

## Usage

``` r
# S3 method for class 'rma.uni'
reliability_profile(x, ...)

# S3 method for class 'rma'
reliability_profile(x, ...)
```

## Arguments

- x:

  An `rma.uni` object, from
  [`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html)
  or
  [`metafor::rma.uni()`](https://wviechtb.github.io/metafor/reference/rma.uni.html).

- ...:

  Must be empty.

## Value

A `reliability_profile` object whose source names the method of the fit.

## Details

The method reads the elements `yi`, `vi`, `tau2`, `method` and `slab` of
the fit. With moderators, \\\tau^2\\ is the residual heterogeneity, the
variance of the true effects around the meta-regression, and the
reliabilities are computed from it. In metafor's notation \\\tau^2\\ is
what this package writes \\\sigma\_\tau^2\\. A fit whose \\\tau^2\\ is
estimated as 0, as often happens with few studies, gives the limit of
the profile at `sigma_tau = 0` and a message; the posterior has no such
boundary.

[`map_tau2_interval()`](https://joonho112.github.io/poolinfo/reference/map_tau2_interval.md)
maps the confidence interval of the fit through the diagnostics. Other
fits of metafor, such as those of
[`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html)
with several variance components, are refused with an error: the
diagnostics need one between-study variance and known sampling
variances. Fixed-, equal- and common-effect fits (`method = "FE"`,
`"EE"` or `"CE"`) are also refused. Refit a random-effects model, or
supply the estimates and sampling variances directly for a new analysis.

## See also

[`map_tau2_interval()`](https://joonho112.github.io/poolinfo/reference/map_tau2_interval.md),
and
[`vignette("meta-analysis")`](https://joonho112.github.io/poolinfo/articles/meta-analysis.md)
for an example.

## Examples

``` r
fit <- metafor::rma(yi = estimate, sei = se, slab = school, data = eight_schools)
reliability_profile(fit)
#> At `sigma_tau = 0` the profile is a limit: every reliability is 0.
#> ℹ The concentration count is its limit as sigma_tau goes to 0.
#> Reliability profile of 8 sites at sigma_tau = 0 (metafor REML estimate)
#> 
#>   informativeness I            0.000
#>   mean reliability             0.000
#>   dispersion H                 0.0000
#>   logit variance V             0.230
#>   concentration count J_eff    6.8 of 8 (limit as sigma_tau goes to 0)
#> 
#> At sigma_tau = 0 every reliability is 0.
reliability_posterior(fit)
#> Reliability profile of 8 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#>                              mean     95% interval
#>   informativeness I          0.225    0.000    0.745
#>   mean reliability           0.230    0.000    0.735
#>   dispersion H               0.0143   0.0000   0.0283
#>   concentration count J_eff  7.2      6.8      7.9      of 8
#>   logit variance V           0.230    (fixed by the standard errors)
#>   sigma_tau                  6.58     0.243    20.8
#>   effect mean                7.93     -2.09    18.2
#> 
#>   Pr(I > 0.5)                0.138
#>   I at E(sigma_tau^2)        0.338
```
