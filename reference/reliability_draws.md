# Diagnostics from posterior draws of the between-site SD

Computes the diagnostics of the reliability profile at each posterior
draw of the between-site standard deviation from another fit, such as a
model fitted with Stan or brms, and summarizes them with their Monte
Carlo error.

## Usage

``` r
reliability_draws(
  x,
  ...,
  se = NULL,
  se2 = NULL,
  site = NULL,
  estimate = NULL,
  data = NULL,
  mean = NULL,
  variable = NULL,
  model_label = NULL,
  scale_label = NULL,
  target_label = NULL,
  level = 0.95,
  threshold = 0.5
)
```

## Arguments

- x:

  Posterior draws of the between-site standard deviation: a numeric
  vector, an object that the posterior package can convert to draws
  (from Stan, CmdStan or brms, among others), or a fit of brms, from
  which the standard errors, the estimates and the sites are taken too.

- ...:

  Must be empty. The arguments after `x` are named, because an unnamed
  vector could hold standard errors or sampling variances.

- se, se2:

  Standard errors or sampling variances of the site estimates, as in
  [`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md).
  Give one of them.

- site:

  Optional labels of the sites.

- estimate:

  Optional estimates of the sites. With `mean`, the site table has the
  posterior means and intervals of the site effects.

- data:

  Optional data frame in which `se`, `se2`, `site` and `estimate` are
  found, as in
  [`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md).

- mean:

  Optional draws of the mean of the effects, in the same order as the
  draws of `x`: a numeric vector, or with a draws object the name of its
  variable.

- variable:

  Name of the variable of the standard deviation in a draws object, a
  string. `NULL`, the default, takes the one variable whose name starts
  with `sd_` or is `sigma_tau`, `tau` or `sigma`; when there is not
  exactly one such variable, `variable` must be given.

- model_label:

  A description of the model and the prior that gave the draws, a
  string, printed with the summaries.

- scale_label, target_label:

  Optional strings that name the scale of the effects and the estimand
  at each site, as in
  [`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md).
  They are printed and used by
  [`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md)
  and
  [`report_table()`](https://joonho112.github.io/poolinfo/reference/report_table.md).

- level:

  Probability content of the posterior intervals, a number between 0 and
  1.

- threshold:

  Value of the informativeness whose posterior probability of being
  exceeded is reported, a number between 0 and 1.

## Value

An object of class `reliability_draws`, a list whose main components are

- `summary`: a data frame with the mean, standard deviation, interval,
  Monte Carlo standard error (`mcse`) and effective sample size of the
  draws (`ess`) of each diagnostic, of `sigma_tau` and `sigma_tau^2`,
  and of the mean when its draws are given.

- `draws`: the draws of `sigma_tau` and of the mean.

- `probability` and `probability_mcse`: the proportion of draws with
  \\I\\ above `threshold`, and its Monte Carlo standard error.

- `plugin`: the diagnostics at the mean of the draws of `sigma_tau^2`.

- `se2`, `site`, `estimate`, `level`, `threshold`, `model_label` and
  `call`: the inputs.

[`profile_summary()`](https://joonho112.github.io/poolinfo/reference/profile_summary.md),
[`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md),
[`prob_above()`](https://joonho112.github.io/poolinfo/reference/prob_above.md),
[`plugin_gap()`](https://joonho112.github.io/poolinfo/reference/plugin_gap.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) take it.

## Details

External-draw summaries are an experimental extension beyond Lee (2026);
their interfaces may change in a later 0.x release.

The diagnostics depend on the draws only through the between-site
standard deviation, so they apply to any posterior algorithm whose draws
target the stated model and prior (Lee, 2026, Appendix D.3). Nothing is
refitted: the model and the prior are those of the fit that gave the
draws, and `model_label` names them for the printed output.

The posterior means are draw averages and come with their Monte Carlo
standard errors. With the posterior package installed, the effective
sample sizes of the draws
([`posterior::ess_mean()`](https://mc-stan.org/posterior/reference/ess_mean.html)),
which account for the autocorrelation of each chain, give the standard
errors; otherwise the draws are treated as independent, with a message.
For the integrals of the normal model with a flat or proper prior on the
SD,
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
has no Monte Carlo error.

A univariate Gaussian identity-link fit of brms of the form
`y | se(sei, sigma = FALSE) ~ 1 + (1 | group)` is read directly.
`sigma = FALSE` may be omitted. The response, known standard errors and
grouping variable must be direct data columns, with one observation per
group. The standard deviation comes from `sd_<group>__Intercept` and the
common mean from `b_Intercept`. Extra residual variance, other families,
additional group-level components, random slopes, distributional
formulas, autocorrelation, smooths, Gaussian processes and latent
measurement-error terms are not supported.

Gaussian models with fixed moderators can supply variance-based
diagnostics and site reliabilities. Site effect estimates are omitted: a
single common-mean draw cannot represent their site-specific population
predictions. Giving `mean` for such a fit is an error.

The site summaries use the draws of the mean rather than the flat-prior
formula: the posterior mean of site \\j\\'s effect averages \\\mu +
S_j(\hat\tau_j - \mu)\\ over the draws, and its interval is a quantile
of the normal mixture over them.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md),
whose integrals have no Monte Carlo error, and
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
for a brms example.

## Examples

``` r
# Draws of sigma_tau from some fit; here simulated for the example.
set.seed(1)
sigma_draws <- abs(rnorm(4000, 8, 4))
rd <- reliability_draws(sigma_draws, se = eight_schools$se,
                        model_label = "an example")
rd
#> Reliability profile of 8 sites (posterior draws from an example)
#> 4000 draws; Monte Carlo standard errors of the means in parentheses.
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#>                              mean                95% interval
#>   informativeness I          0.303    (0.0029)   0.005    0.637
#>   mean reliability           0.309    (0.0029)   0.005    0.631
#>   dispersion H               0.0201   (0.00014)  0.0005   0.0284
#>   concentration count J_eff  7.3      (0.0043)   6.8      7.8      of 8
#>   logit variance V           0.230               (fixed by the standard errors)
#>   sigma_tau                  8.10     (0.063)    0.844    16.1
#> 
#>   Pr(I > 0.5)                0.162 (0.0057)
#>   I at E(sigma_tau^2)        0.355
```
