# Reliability profile at a fixed between-site SD

Computes the site reliabilities and the diagnostics of their profile at
a given between-site standard deviation: a planning value, the true
value of a simulation, or an estimate from another fit.

## Usage

``` r
reliability_profile(x = NULL, ...)

# Default S3 method
reliability_profile(
  x = NULL,
  ...,
  se = NULL,
  se2 = NULL,
  sigma_tau,
  site = NULL,
  estimate = NULL,
  data = NULL,
  scale_label = NULL,
  target_label = NULL
)

# S3 method for class 'data.frame'
reliability_profile(
  x,
  ...,
  se = NULL,
  se2 = NULL,
  sigma_tau,
  site = NULL,
  estimate = NULL,
  scale_label = NULL,
  target_label = NULL
)

# S3 method for class 'reliability_profile'
update(object, ..., sigma_tau = NULL, scale_label = NULL, target_label = NULL)

# S3 method for class 'reliability_posterior'
reliability_profile(x, ...)

# S3 method for class 'reliability_draws'
reliability_profile(x, ...)

# S3 method for class 'reliability_profile'
format(x, ..., use_cli = NULL)

# S3 method for class 'reliability_profile'
print(x, ..., use_cli = NULL)

# S3 method for class 'reliability_profile'
summary(object, ...)

# S3 method for class 'summary.reliability_profile'
format(x, ..., use_cli = NULL)

# S3 method for class 'summary.reliability_profile'
print(x, ..., use_cli = NULL)

# S3 method for class 'reliability_profile'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  Usually left out, with the standard errors given as `se` or `se2`. It
  can also be a data frame, in which `se`, `se2`, `site` and `estimate`
  name columns; a posterior fit or posterior draws, for the profile at
  the posterior mean of \\\sigma\_\tau^2\\; or a fit of metafor, ebnm or
  multisiteDGP (see
  [metafor-methods](https://joonho112.github.io/poolinfo/reference/metafor-methods.md),
  [ebnm-methods](https://joonho112.github.io/poolinfo/reference/ebnm-methods.md)
  and
  [multisitedgp-methods](https://joonho112.github.io/poolinfo/reference/multisitedgp-methods.md)).

- ...:

  Must be empty.

- se:

  Standard errors of the site estimates, a numeric vector of positive
  numbers. With a data frame `x`, an unquoted column name. Give `se` or
  `se2`, not both.

- se2:

  Sampling variances of the site estimates, as `se`.

- sigma_tau:

  Standard deviation of the true effects across sites, a number of at
  least 0.

- site:

  Optional labels of the sites, one for each site and all distinct. With
  a data frame `x`, an unquoted column name.

- estimate:

  Optional estimates of the sites, a numeric vector. When they are
  given, the site table also has the posterior mean and standard
  deviation of each site's effect at `sigma_tau`.

- data:

  Optional data frame in which `se`, `se2`, `site` and `estimate` are
  found, as in
  [`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md);
  the same as giving the data frame as `x`.

- scale_label, target_label:

  Optional strings that name the scale of the effects, such as
  `"risk difference"` (read as "on the risk difference scale"), and the
  estimand at each site, such as
  `"the block-adjusted intention-to-treat contrast for course passing"`.
  They appear in the printed output and in
  [`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md).

- object:

  A profile made by `reliability_profile()`.

- use_cli:

  Whether to print the headings in bold with the cli package, when it is
  installed. `NULL`, the default, takes
  `getOption("poolinfo.use_cli", TRUE)`.

- row.names, optional:

  Ignored.

## Value

An object of class `reliability_profile`, a list with components

- `se2`, `site` and `estimate` (or `NULL`): the data.

- `X`: the mean-model matrix inherited from a fitted model, or `NULL`
  for a common mean.

- `sigma_tau` and `source`: the between-site SD and where it came from.

- `S`: the reliabilities.

- `stats`: a named vector of the diagnostics: `I`, `mean`, `harmonic`,
  `H`, `V`, `rho3` (the third absolute moment of the centered logits),
  `J_eff`, `gap` (the mean reliability minus \\I\\), `benchmark_risk`
  (see
  [`benchmark_risk()`](https://joonho112.github.io/poolinfo/reference/benchmark_risk.md)),
  `trace` (the sum of the reliabilities), `v_typ` (the typical variance
  of
  [`typical_variance()`](https://joonho112.github.io/poolinfo/reference/typical_variance.md))
  and `R_typ` (the reliability at that variance).

- `local`: the table of
  [`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md).

- `labels`: the scale and estimand labels.

- `call` and `version`: the call and the version of the package that
  made the profile.

It has methods for [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (the site
table), [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`update()`](https://rdrr.io/r/stats/update.html).

## Details

With sampling variances \\v_j\\ and the between-site standard deviation
\\\sigma\_\tau\\, the reliability of site \\j\\ is \\S_j =
\sigma\_\tau^2/(\sigma\_\tau^2 + v_j)\\. The profile's diagnostics are
the informativeness \\I\\
([`kl_centroid()`](https://joonho112.github.io/poolinfo/reference/kl_centroid.md)),
the mean reliability, the harmonic centroid, the dispersion \\H\\
([`kl_dispersion()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md)),
the logit variance \\V\\
([`logit_moments()`](https://joonho112.github.io/poolinfo/reference/logit_moments.md)),
the concentration count \\J\_{\mathrm{eff}}\\
([`concentration_count()`](https://joonho112.github.io/poolinfo/reference/concentration_count.md))
and the terms of
[`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md).
[`summary()`](https://rdrr.io/r/base/summary.html) adds the benchmark
risk, the sum of the reliabilities and the typical variance of
[`typical_variance()`](https://joonho112.github.io/poolinfo/reference/typical_variance.md).

This is the profile at one value of `sigma_tau`. When `sigma_tau` is
estimated from the same data,
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
averages the diagnostics over its posterior distribution instead, which
is what we recommend for reporting. Profiles obtained from a posterior
or a metafor fit with moderators retain its mean-model matrix. Their
site effects pool toward the site's fitted mean, with uncertainty in the
mean-model coefficients integrated out under a flat prior.

At `sigma_tau = 0` every reliability is 0 and the profile is a limit:
\\I\\, the mean reliability and \\H\\ are 0, and the concentration count
is its limit \\(\sum_j v_j^{-1})^2 / \sum_j v_j^{-2}\\, and a message is
shown.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`profile_summary()`](https://joonho112.github.io/poolinfo/reference/profile_summary.md)
and
[`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md)
for the numbers as data frames,
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
for the posterior over `sigma_tau`, and
[`vignette("reliability-profiles")`](https://joonho112.github.io/poolinfo/articles/reliability-profiles.md)
for what the diagnostics describe.

## Examples

``` r
se <- c(0.10, 0.12, 0.15, 0.20, 0.25, 0.30)
prof <- reliability_profile(se = se, sigma_tau = 0.2)
prof
#> Reliability profile of 6 sites at sigma_tau = 0.2 (supplied)
#> 
#>   informativeness I            0.571
#>   mean reliability             0.562
#>   dispersion H                 0.0720
#>   logit variance V             0.609
#>   concentration count J_eff    5.5 of 6
#> 
#> Site reliabilities range from 0.31 to 0.80.
summary(prof)
#> Reliability profile of 6 sites at sigma_tau = 0.2 (supplied)
#> 
#> Centroids
#>   harmonic centroid            0.501    (harmonic mean)
#>   informativeness I            0.571    (logit mean)
#>   mean reliability             0.562    (arithmetic mean)
#> 
#> Dispersion and concentration
#>   dispersion H                 0.0720
#>   logit variance V             0.609
#>   concentration count J_eff    5.5 of 6
#> 
#> Second-order terms (Lee, 2026)
#>                                                           leading term
#>                           exact      leading    bound     decides the sign
#>   H                       0.0720     0.0746     0.00933   -
#>   mean reliability - I   -0.00922   -0.0107     0.0121    no
#> 
#> Related summaries
#>   benchmark risk sigma_tau^2 (1 - mean S_j)    0.0175
#>   sum of reliabilities                         3.4
#>   typical variance (Higgins-Thompson)          0.0251   R_typ 0.614

# The same from a data frame, with the columns named
dat <- data.frame(school = letters[1:6], y = c(0.3, 0.1, 0.25, -0.05, 0.4, 0),
                  s = se)
reliability_profile(dat, se = s, estimate = y, site = school, sigma_tau = 0.2)
#> Reliability profile of 6 sites at sigma_tau = 0.2 (supplied)
#> 
#>   informativeness I            0.571
#>   mean reliability             0.562
#>   dispersion H                 0.0720
#>   logit variance V             0.609
#>   concentration count J_eff    5.5 of 6
#> 
#> Site reliabilities range from 0.31 to 0.80.

# The profile at another sigma_tau
update(prof, sigma_tau = 0.1)
#> Reliability profile of 6 sites at sigma_tau = 0.1 (supplied)
#> 
#>   informativeness I            0.250
#>   mean reliability             0.276
#>   dispersion H                 0.0565
#>   logit variance V             0.609
#>   concentration count J_eff    4.7 of 6
#> 
#> Site reliabilities range from 0.10 to 0.50.
```
