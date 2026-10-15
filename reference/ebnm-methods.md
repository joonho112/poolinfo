# Profiles of ebnm fits with a normal prior

[`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
of a fit of the ebnm package with a normal prior gives the profile at
the standard deviation of the fitted prior.

## Usage

``` r
# S3 method for class 'ebnm'
reliability_profile(x, ...)
```

## Arguments

- x:

  An `ebnm` object fitted with `prior_family = "normal"`.

- ...:

  Must be empty.

## Value

A `reliability_profile` object whose source is the ebnm fit.

## Details

The method reads the elements `data$x`, `data$s` and `fitted_g` of the
fit, where `fitted_g` must be a normal distribution (a `normalmix` with
one component). ebnm's normal prior has its mean fixed at 0 unless the
fit used `mode = "estimate"`; the reliabilities depend only on the
standard deviation, but the site table's posterior means use the
precision-weighted mean and so differ from ebnm's own when the mode was
fixed. A prior of another family is not a normal model and is refused.

## Examples

``` r
fit <- ebnm::ebnm(eight_schools$estimate, eight_schools$se, prior_family = "normal",
                  mode = "estimate")
reliability_profile(fit)
#> Reliability profile of 8 sites at sigma_tau = 0.00244 (ebnm normal prior)
#> 
#>   informativeness I            0.000
#>   mean reliability             0.000
#>   dispersion H                 0.0000
#>   logit variance V             0.230
#>   concentration count J_eff    6.8 of 8
#> 
#> Site reliabilities range from 0.00 to 0.00.
```
