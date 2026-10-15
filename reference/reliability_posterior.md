# Posterior summaries of the reliability profile

Computes the posterior distribution of the informativeness, the mean
reliability, the dispersion and the concentration count of a set of site
estimates with standard errors, under the normal hierarchical model, and
the posterior means and intervals of the site effects.

## Usage

``` r
reliability_posterior(
  estimate,
  ...,
  se = NULL,
  se2 = NULL,
  site = NULL,
  data = NULL,
  mean_model = ~1,
  prior = sd_prior("flat"),
  level = 0.95,
  threshold = 0.5,
  scale_label = NULL,
  target_label = NULL,
  na_action = c("fail", "omit"),
  control = posterior_control()
)

# S3 method for class 'reliability_posterior'
format(x, ..., use_cli = NULL)

# S3 method for class 'reliability_posterior'
print(x, ..., use_cli = NULL)

# S3 method for class 'reliability_posterior'
summary(object, ...)

# S3 method for class 'summary.reliability_posterior'
format(x, ..., use_cli = NULL)

# S3 method for class 'summary.reliability_posterior'
print(x, ..., use_cli = NULL)

# S3 method for class 'reliability_posterior'
update(
  object,
  ...,
  prior = NULL,
  level = NULL,
  threshold = NULL,
  control = NULL
)

# S3 method for class 'reliability_posterior'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- estimate:

  The site estimates: a numeric vector, or a column of `data`. It can
  also be a random-effects fit of metafor (`rma.uni`) or a data set
  simulated by multisiteDGP, which carry the estimates and their
  sampling variances (see
  [metafor-methods](https://joonho112.github.io/poolinfo/reference/metafor-methods.md)
  and
  [multisitedgp-methods](https://joonho112.github.io/poolinfo/reference/multisitedgp-methods.md)).

- ...:

  Must be empty. The arguments after `estimate` are named, because an
  unnamed vector could hold standard errors or sampling variances.

- se, se2:

  The standard errors of the estimates, or their sampling variances (the
  squared standard errors): positive numbers, or a column of `data`.
  Give one of the two.

- site:

  Optional labels of the sites, all distinct, used in the tables and
  plots.

- data:

  A data frame in which `estimate`, `se`, `se2`, `site` and the
  variables of `mean_model` are found. A column takes precedence over a
  variable of the same name in the calling environment; write `.env$x`
  for the variable.

- mean_model:

  A one-sided formula for the mean of the effects, as in a
  meta-regression; `~ 1`, the default, is a common mean.

- prior:

  The prior of the between-site standard deviation, made by
  [`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md);
  flat by default. In [`update()`](https://rdrr.io/r/stats/update.html),
  `NULL` keeps the prior of the fit, as it does for `level`, `threshold`
  and `control`.

- level:

  Probability content of the posterior intervals, a number between 0 and
  1.

- threshold:

  The value of the informativeness whose posterior probability of being
  exceeded is reported, between 0 and 1. See
  [`prob_above()`](https://joonho112.github.io/poolinfo/reference/prob_above.md).

- scale_label, target_label:

  Optional strings that name the scale of the effects, such as
  `"risk difference"` (read as "on the risk difference scale"), and the
  estimand at each site, such as
  `"the block-adjusted intention-to-treat contrast for course passing"`.
  They appear in the printed output and in
  [`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md).

- na_action:

  `"fail"` (the default) stops at a site with a missing value; `"omit"`
  leaves such sites out, with a message that names them.

- control:

  Settings of the numerical integration, from
  [`posterior_control()`](https://joonho112.github.io/poolinfo/reference/posterior_control.md).

- x, object:

  A fit of `reliability_posterior()`.

- use_cli:

  Whether to print the headings in bold with the cli package, when it is
  installed. `NULL`, the default, takes
  `getOption("poolinfo.use_cli", TRUE)`.

- row.names, optional:

  Ignored.

## Value

`reliability_posterior()` and
[`update()`](https://rdrr.io/r/stats/update.html) return an object of
class `reliability_posterior`, a list whose main components are

- `summary`: a data frame with the posterior mean, standard deviation
  and interval of each diagnostic, of `sigma_tau`, `sigma_tau^2` and of
  the mean or the coefficients.

- `quantiles`: the lower end, median and upper end of the same
  quantities.

- `probability`: the posterior probability that \\I\\ exceeds
  `threshold`.

- `plugin`: the diagnostics at the posterior mean of `sigma_tau^2`.

- `sites`: the posterior means of the site reliabilities and site
  effects.

- `checks`: the checks of the numerical integration.

- `data`, `prior`, `control`, `level`, `threshold`, `labels` and `call`:
  the inputs.

The accessors
[`profile_summary()`](https://joonho112.github.io/poolinfo/reference/profile_summary.md),
[`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md),
[`prob_above()`](https://joonho112.github.io/poolinfo/reference/prob_above.md),
[`plugin_gap()`](https://joonho112.github.io/poolinfo/reference/plugin_gap.md),
[`functional_density()`](https://joonho112.github.io/poolinfo/reference/functional_density.md),
[`profile_draws()`](https://joonho112.github.io/poolinfo/reference/profile_draws.md)
and
[`check_numerics()`](https://joonho112.github.io/poolinfo/reference/check_numerics.md)
return these numbers as data frames and vectors, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows
[`plot_dashboard()`](https://joonho112.github.io/poolinfo/reference/plot_dashboard.md).
[`summary()`](https://rdrr.io/r/base/summary.html) returns an object
that prints the longer summary, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the site table.

## Details

The working model is \$\$\hat\tau_j \mid \tau_j \sim N(\tau_j, v_j),
\qquad \tau_j \sim N(x_j^\top\beta, \sigma\_\tau^2),\$\$ independently
over sites, with the sampling variances \\v_j\\ taken as known. The
coefficients \\\beta\\ (by default the common mean) have a flat prior
and are integrated out, which leaves a posterior for `sigma_tau` alone
(Lee, 2026, Appendix D). The default prior of `sigma_tau` is flat too,
as in Lee (2026). Each summary is a one-dimensional integral over this
posterior, computed numerically to about 10 digits, so there is no Monte
Carlo error.

The diagnostics are functions of `sigma_tau` and the \\v_j\\. The
informativeness \\I\\, the mean reliability, the harmonic centroid and
the concentration count increase with `sigma_tau`. The posterior
intervals are central, with probability \\(1 - \mathrm{level})/2\\ below
and above them, so for these four quantities they are the images of the
interval of `sigma_tau`. The dispersion \\H\\ first rises and then
falls, and its quantiles take both branches into account. The logit
variance \\V\\ depends on the \\v_j\\ only and has no posterior
distribution.

The summaries are conditional on the supplied standard errors, which are
treated as exact, and every site counts the same in them, whatever its
size. The last two lines of the printed header state both.

Under the flat prior the posterior is proper with at least \\p + 2\\
sites, where \\p\\ is the number of columns of the mean model (3 sites
for a common mean). The posterior mean of `sigma_tau` needs at least
\\p + 3\\ sites and that of `sigma_tau^2` at least \\p + 4\\. A
posterior standard deviation needs the next moment as well: at least
\\p + 4\\ sites for `sigma_tau` and the mean, and \\p + 6\\ for
`sigma_tau^2`. The posterior mean of the common mean or a mean-model
coefficient requires \\E(\sigma\_\tau)\\ to exist, hence at least \\p +
3\\ sites under the flat prior. With fewer sites its mean is not
defined, although a central posterior interval is available. Undefined
moments are `NA`, and the printed fit marks them as not defined.
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md)
gives the conditions under the other priors. The paper's results are for
a common mean and the flat prior. Mean models with covariates and proper
SD priors are experimental extensions beyond Lee (2026); their
interfaces may change in a later 0.x release.
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
describes their computation.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
for a profile at a fixed `sigma_tau`;
[`reliability_draws()`](https://joonho112.github.io/poolinfo/reference/reliability_draws.md)
for posterior draws from other software;
[`compare_profiles()`](https://joonho112.github.io/poolinfo/reference/compare_profiles.md)
and
[`prior_sensitivity()`](https://joonho112.github.io/poolinfo/reference/prior_sensitivity.md)
for several analyses;
[`vignette("poolinfo")`](https://joonho112.github.io/poolinfo/articles/poolinfo.md)
for a worked analysis and
[`vignette("posterior-summaries")`](https://joonho112.github.io/poolinfo/articles/posterior-summaries.md)
for the choices behind the summaries.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
fit
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
profile_summary(fit)
#>      quantity        mean          sd         lower        upper
#> 1           I  0.22529451   0.2173025  3.997828e-04   0.74517275
#> 2        mean  0.22972844   0.2150459  4.438244e-04   0.73510308
#> 3    harmonic  0.21033305   0.2084266  3.547019e-04   0.72178667
#> 4           H  0.01428773   0.0100617  4.407725e-05   0.02832835
#> 5       J_eff  7.17004472   0.3177007  6.818171e+00   7.86891743
#> 6           V  0.23016500   0.0000000  2.301650e-01   0.23016500
#> 7   sigma_tau  6.57548335   5.6503793  2.426961e-01  20.75245613
#> 8  sigma_tau2 75.16376789 162.2015906  5.890141e-02 430.66443537
#> 9 effect_mean  7.93237502   5.1783948 -2.090220e+00  18.21814307
update(fit, prior = sd_prior("half_normal", scale = 10))
#> Reliability profile of 8 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: half-normal on sigma_tau with scale 10, flat on the mean
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#>                              mean     95% interval
#>   informativeness I          0.157    0.000    0.570
#>   mean reliability           0.162    0.000    0.568
#>   dispersion H               0.0118   0.0000   0.0282
#>   concentration count J_eff  7.1      6.8      7.7      of 8
#>   logit variance V           0.230    (fixed by the standard errors)
#>   sigma_tau                  4.84     0.189    14.0
#>   effect mean                7.86     -1.27    17.1
#> 
#>   Pr(I > 0.5)                0.049
#>   I at E(sigma_tau^2)        0.203
```
