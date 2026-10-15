# Posterior summaries of the diagnostics

[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
reports a posterior mean and a 95% posterior interval for each
diagnostic, and the posterior probability that the informativeness $`I`$
exceeds 1/2. In this vignette we explain where these numbers come from
and which choices lie behind them. We start with the model and its
prior, compare a posterior mean with a plug-in value, and explain the
interval for the dispersion $`H`$. We then come back to the prior for
trials with few sites, and end with posterior draws, for summaries that
the package does not compute. Our examples are the eight schools of
Rubin ([1981](#ref-rubin1981)) and the passing outcome of the Online
Credit Recovery Study ([Rickles et al. 2024](#ref-rickles2024)).

``` r

library(poolinfo)
```

## The model and the prior

We work with the normal hierarchical model in which site $`j`$’s
estimate $`\hat\tau_j`$ has a known sampling variance $`v_j`$,

``` math
\hat\tau_j \mid \tau_j \sim N(\tau_j, v_j), \qquad
\tau_j \sim N(\mu, \sigma_\tau^2), \qquad j = 1, \dots, J,
```

independently over sites. By default the prior is flat on the mean
$`\mu`$ and on the between-site standard deviation $`\sigma_\tau`$. A
prior flat on $`\sigma_\tau`$ is not flat on the variance. As a density
in $`\sigma_\tau^2`$ it is proportional to $`1/\sigma_\tau`$, so it puts
less weight on large variances than a flat prior on $`\sigma_\tau^2`$
would ([Gelman 2006](#ref-gelman2006)). Integrating $`\mu`$ out leaves a
posterior for $`\sigma_\tau`$ alone, known up to a constant, and
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
gives its formula. Every diagnostic depends only on $`\sigma_\tau`$ and
the $`v_j`$, so its posterior mean, its quantiles and the probability
that it exceeds a value are one-dimensional integrals over that
posterior. Because the package computes them by quadrature rather than
simulation, they have no Monte Carlo error.

Our first example, `eight_schools`, has the estimated effects of
coaching on the SAT-V in eight schools, with their standard errors. We
fit the posterior.

``` r

schools <- reliability_posterior(estimate, se = se, site = school,
                                 data = eight_schools,
                                 scale_label = "SAT-V point")
schools
#> Reliability profile of 8 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Scale of the effects: SAT-V point
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

The posterior of $`\sigma_\tau`$ is wide. Its 95% interval runs from
0.24 to 20.8 SAT-V points, and the interval for $`I`$ runs from 0.00 to
0.75.

In the last block, Pr(I \> 0.5) is the posterior probability that
$`\sigma_\tau^2`$ exceeds $`\mathrm{GM}(v)`$, the geometric mean of the
sampling variances. The threshold 1/2 is where the leading term of
$`\bar S - I`$ changes sign ([Lee 2026](#ref-lee2026), Proposition 2),
as
[`vignette("poolinfo")`](https://joonho112.github.io/poolinfo/articles/poolinfo.md)
explains.

[`plot_sigma_scale()`](https://joonho112.github.io/poolinfo/reference/plot_sigma_scale.md)
plots the posterior density of $`\sigma_\tau`$ with a second axis for
$`I`$, and shades the region where $`I`$ exceeds 1/2.

``` r

plot_sigma_scale(schools)
```

![The posterior density of sigma_tau for the eight schools, highest at 0
and falling, with a second axis for the informativeness and the region
where it exceeds 1/2
shaded.](posterior-summaries_files/figure-html/sigma-scale-1.png)

The density is highest at 0 and falls steadily. The shading starts at
$`\sigma_\tau = \sqrt{\mathrm{GM}(v)} = 12.1`$ points, where
$`I = 1/2`$, and 0.138 of the probability lies beyond it. Because $`I`$
increases with $`\sigma_\tau`$, each quantile of $`I`$ is $`I`$ at the
same quantile of $`\sigma_\tau`$, which the top axis lets you read off.

## Posterior means and plug-in values

The printed fit gives the posterior mean of $`I`$ in its first row and,
in its last line, $`I`$ at the posterior mean of $`\sigma_\tau^2`$. That
is the value you would get by putting this estimate of $`\sigma_\tau^2`$
into the formula for $`I`$.
[`plugin_gap()`](https://joonho112.github.io/poolinfo/reference/plugin_gap.md)
gives both, and their difference, for four diagnostics.

``` r

print(plugin_gap(schools), digits = 3)
#>   functional posterior_mean at_mean_sigma2 difference
#> 1          I         0.2253         0.3379     0.1126
#> 2       mean         0.2297         0.3464     0.1167
#> 3          H         0.0143         0.0249     0.0106
#> 4      J_eff         7.1700         7.3631     0.1930
```

For the eight schools the two values of $`I`$ are far apart, 0.225
against 0.338. The posterior mean averages $`I`$ over the values of
$`\sigma_\tau`$ that the data allow. The plug-in value first averages
$`\sigma_\tau^2`$ and then applies $`I`$, and the long right tail pulls
the posterior mean of $`\sigma_\tau^2`$, 75.2, far above its median,
27.4. Whatever the data, $`I`$ is a concave function of
$`\sigma_\tau^2`$. By Jensen’s inequality, then, its posterior mean is
at most its value at the posterior mean of $`\sigma_\tau^2`$, with
equality only when the posterior of $`\sigma_\tau^2`$ is a point mass
([Lee 2026](#ref-lee2026), Proposition 3). The same holds for every
$`S_j`$, and so for $`\bar S`$.

Proposition 3 gives no such order for $`H`$ and $`J_{\mathrm{eff}}`$,
and the positive differences in their rows are a property of these data.
For the passing outcome of the credit recovery study, with 32 sites, the
difference for $`H`$ is negative.

``` r

passing <- reliability_posterior(estimate, se = se, data = ocrs)
print(plugin_gap(passing), digits = 3)
#>   functional posterior_mean at_mean_sigma2 difference
#> 1          I         0.6709         0.6911   2.02e-02
#> 2       mean         0.6635         0.6823   1.88e-02
#> 3          H         0.0253         0.0253  -2.04e-05
#> 4      J_eff        31.1704        31.2816   1.11e-01
```

Here $`I`$ and $`\bar S`$ again come out larger at the posterior mean of
$`\sigma_\tau^2`$, as they must, while $`H`$ comes out smaller, by
2.0e-05.

A third number is common in practice, $`I`$ at a point estimate of
$`\sigma_\tau^2`$ such as the REML estimate. For the eight schools that
estimate is 0, which makes every reliability 0.
[`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
computes the profile at any value of $`\sigma_\tau`$, or at the estimate
of a metafor fit, but such a profile carries no uncertainty about
$`\sigma_\tau`$.

## Intervals

The posterior intervals are central, with 2.5% of the posterior
probability below the lower end and 2.5% above the upper end. We use
central intervals because they carry over through increasing functions.
$`I`$, $`\bar S`$ and $`J_{\mathrm{eff}}`$ increase with
$`\sigma_\tau`$, so the ends of their intervals are their values at the
2.5% and 97.5% quantiles of $`\sigma_\tau`$. A highest-density interval
does not carry over in general, since a nonlinear change of scale moves
the region of highest density. For the eight schools it would also run
from 0, where the density of $`\sigma_\tau`$ is highest, to the 95%
quantile.

The dispersion $`H`$ needs more care. It is 0 at $`\sigma_\tau = 0`$,
peaks where $`\bar S = I`$, and returns toward 0 as $`\sigma_\tau`$
grows, so each value of $`H`$ below the peak comes from two values of
$`\sigma_\tau`$. The package adds the posterior probabilities of the two
branches to get the distribution function of $`H`$, and solves it for
the quantiles. For the eight schools the interval for $`H`$ runs from
0.00004 to 0.02833, just below the largest value that $`H`$ can take
with these standard errors, 0.02837. The upper end is so close to the
largest value because $`H`$ is flat near its peak. Every $`\sigma_\tau`$
from 12.1 to 13.1 points gives a value of $`H`$ above the upper end, and
these values of $`\sigma_\tau`$ carry the remaining 2.5% of the
probability.

[`functional_density()`](https://joonho112.github.io/poolinfo/reference/functional_density.md)
returns the posterior density of a diagnostic on a grid, by a change of
variables from the density of $`\sigma_\tau`$, and
[`plot_functionals()`](https://joonho112.github.io/poolinfo/reference/plot_functionals.md)
plots it with the posterior mean, the plug-in value and the 95%
interval.

``` r

plot_functionals(schools)
```

![The posterior density of the informativeness for the eight schools,
unbounded near 0, with the posterior mean as a solid line, the value at
the posterior mean of sigma_tau squared as a dashed line and the 95%
interval shaded.](posterior-summaries_files/figure-html/density-1.png)

The solid line marks the posterior mean of $`I`$, 0.225, the dashed line
the plug-in value, 0.338, and the shading the interval up to 0.75. Near
$`I = 0`$ the density is unbounded, and the plot cuts its vertical axis.
Under the flat prior, as under every prior that
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md)
offers, the posterior density of $`\sigma_\tau`$ is positive at 0. Close
to 0 the informativeness is nearly $`\sigma_\tau^2/\mathrm{GM}(v)`$, so
the density of $`I`$ grows like $`I^{-1/2}`$ as $`I`$ approaches 0. A
gray line marks the threshold 1/2.

[`prob_above()`](https://joonho112.github.io/poolinfo/reference/prob_above.md)
computes the posterior probability that $`I`$, $`\bar S`$,
$`J_{\mathrm{eff}}`$ or $`\sigma_\tau`$ itself exceeds a value. We ask
for the probability that $`\sigma_\tau`$ exceeds 5 SAT-V points.

``` r

prob_above(schools, threshold = 5, which = "sigma_tau")
#> [1] 0.5194775
```

It is 0.52, a little over a half.

## Few sites and the prior

With few sites the posterior of $`\sigma_\tau`$ is wide, and the prior
matters. Under the flat prior, the package’s default and the prior of
Lee ([2026](#ref-lee2026)), the posterior is proper only with at least 3
sites, the posterior mean of $`\sigma_\tau`$ exists only with at least
4, and that of $`\sigma_\tau^2`$ only with at least 5.
[`?sd_prior`](https://joonho112.github.io/poolinfo/reference/sd_prior.md)
gives these counts for a mean model and for the other priors. When a
posterior mean does not exist, the printed fit reports it as not
defined. We fit the first 4 of the eight schools.

``` r

four <- reliability_posterior(estimate, se = se, data = eight_schools[1:4, ])
four
#> Reliability profile of 4 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#>                              mean     95% interval
#>   informativeness I          0.439    0.001    0.976
#>   mean reliability           0.439    0.001    0.975
#>   dispersion H               0.0112   0.0001   0.0197
#>   concentration count J_eff  3.8      3.5      4.0      of 4
#>   logit variance V           0.159    (fixed by the standard errors)
#>   sigma_tau                  18.7     0.465    82.0
#>   effect mean                9.71     -18.6    38.5
#> 
#>   Pr(I > 0.5)                0.437
#>   I at E(sigma_tau^2)        not defined: E(sigma_tau^2) does not exist with 4 sites
```

The posterior mean of $`I`$ is now 0.44, against 0.23 with all 8
schools, and the 95% interval of $`\sigma_\tau`$ reaches 82 points. For
both sets the REML estimate of $`\sigma_\tau`$ is 0, but four schools
say less about $`\sigma_\tau`$ than eight. Under the flat prior the
posterior density of $`\sigma_\tau`$ falls off like
$`\sigma_\tau^{1-J}`$, slowly when $`J = 4`$, and the posterior
probability that $`\sigma_\tau`$ exceeds 20 points is 0.27 with 4
schools against 0.03 with 8. The posterior mean of $`\sigma_\tau^2`$
does not exist with 4 sites, so the printout has no plug-in value. Under
a half-normal prior with scale 10 the two sets differ as well, with
posterior means of $`I`$ of 0.20 and 0.16.

Proper priors on $`\sigma_\tau`$ are a common alternative when there are
few sites, such as half-normal or half-t priors with a scale on the
scale of the effects ([Gelman 2006](#ref-gelman2006)).
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md)
makes them, and
[`prior_sensitivity()`](https://joonho112.github.io/poolinfo/reference/prior_sensitivity.md)
refits a posterior under each and collects the fits:

``` r

prior_sensitivity(schools, list(
  `half-normal, scale 10` = sd_prior("half_normal", scale = 10),
  `half-t, scale 10, 3 df` = sd_prior("half_t", scale = 10, df = 3),
  `half-Cauchy, scale 5` = sd_prior("half_cauchy", scale = 5)
))
#> Comparison of 4 analyses
#> Intervals are central posterior intervals.
#> Each analysis is conditional on its own standard errors.
#> 
#>   analysis                J  informativeness I     mean reliability
#>   flat                    8  0.225 [0.000, 0.745]  0.230
#>   half-normal, scale 10   8  0.157 [0.000, 0.570]  0.162
#>   half-t, scale 10, 3 df  8  0.156 [0.000, 0.587]  0.161
#>   half-Cauchy, scale 5    8  0.105 [0.000, 0.497]  0.110
#> 
#>   analysis                dispersion H             concentration count J_eff
#>   flat                    0.0143 [0.0000, 0.0283]  7.2
#>   half-normal, scale 10   0.0118 [0.0000, 0.0282]  7.1
#>   half-t, scale 10, 3 df  0.0116 [0.0000, 0.0282]  7.1
#>   half-Cauchy, scale 5    0.0084 [0.0000, 0.0278]  7.0
```

All three put more weight on small values of $`\sigma_\tau`$ than the
flat prior does. They pull the posterior mean of $`I`$ down from 0.23 to
between 0.11 and 0.16, and the interval’s upper end from 0.75 to between
0.50 and 0.59. We would name the prior in any report on a trial this
small, and show the fits under more than one.

A given prior matters less as the number of sites grows. To see this
with the trial and the scale of the effects held fixed, we refit the
passing outcome in its 13 Algebra 1 sites, its 19 English 9 sites and
all 32 sites, each under two half-normal priors. A scale of 0.1 is about
0.8 times the geometric-mean standard error $`\sqrt{\mathrm{GM}(v)}`$ in
each set, as the scale 10 is for the eight schools. A scale of 0.5, a
between-site standard deviation of 50 percentage points, is far beyond
any plausible spread of effects on passing, so that prior is weak. For
each set we compute how far each prior moves the posterior mean of $`I`$
from its value under the flat prior.

``` r

half_normal <- list(`scale 0.1` = sd_prior("half_normal", scale = 0.1),
                    `scale 0.5` = sd_prior("half_normal", scale = 0.5))
subsets <- list(`Algebra 1` = ocrs[ocrs$subject == "Algebra 1", ],
                `English 9` = ocrs[ocrs$subject == "English 9", ],
                all = ocrs)
shift <- sapply(subsets, function(d) {
  fit <- reliability_posterior(estimate, se = se, data = d)
  I <- as.data.frame(prior_sensitivity(fit, half_normal))$I
  setNames(I[-1] - I[1], names(half_normal))
})
data.frame(sites = sapply(subsets, nrow), round(t(shift), 3),
           check.names = FALSE)
#>           sites scale 0.1 scale 0.5
#> Algebra 1    13    -0.118    -0.007
#> English 9    19    -0.072    -0.004
#> all          32    -0.045    -0.002
```

Under the prior with scale 0.1 the posterior mean of $`I`$ falls by
0.118 with 13 sites, 0.072 with 19 and 0.045 with 32. Under the weak
prior it falls by 0.007 or less. So a weak prior hardly matters with 32
sites, but a prior with a scale near the standard errors still moves
$`I`$ by about half of its posterior standard deviation, 0.081.

## Draws

Some summaries need posterior draws, for example those that involve
several site effects at once.
[`profile_draws()`](https://joonho112.github.io/poolinfo/reference/profile_draws.md)
draws $`\sigma_\tau`$ by inverting its posterior distribution function,
then the mean given $`\sigma_\tau`$, and, with `sites = TRUE`, the site
effects given both. With `seed` the draws are reproducible, and the
session’s random number generator is left as it was. We use draws to
find the posterior probability that school A, whose estimate is the
largest, also has the largest effect:

``` r

draws <- profile_draws(schools, ndraws = 4000, seed = 1, sites = TRUE)
effects <- draws[, paste0("effect_", eight_schools$school)]
top_A <- mean(apply(effects, 1, which.max) == 1)
round(c(probability = top_A,
        mcse = sqrt(top_A * (1 - top_A) / nrow(draws))), 4)
#> probability        mcse 
#>      0.2585      0.0069
```

School A has the largest effect in 26% of the draws, more often than any
other school but far from always, because the model pulls the eight
effects toward one another. The draws are independent, so the Monte
Carlo standard error of this proportion is the binomial one, 0.007. Four
times as many draws would halve it.

Draws from other software go to
[`reliability_draws()`](https://joonho112.github.io/poolinfo/reference/reliability_draws.md),
which computes the diagnostics at each draw of $`\sigma_\tau`$ and
reports their Monte Carlo standard errors. It reads a brms fit of this
model directly, and
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
has an example.

## Integration checks

The numerical integration is checked as each fit is made, and
[`check_numerics()`](https://joonho112.github.io/poolinfo/reference/check_numerics.md)
shows the results. If a check falls outside its limit,
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
warns as it builds the fit.
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
explains the checks and the other method of integration,
`posterior_control(engine = "adaptive")`.

## References

Gelman, Andrew. 2006. “Prior Distributions for Variance Parameters in
Hierarchical Models (Comment on Article by Browne and Draper).”
*Bayesian Analysis* 1 (3): 515–34. <https://doi.org/10.1214/06-BA117A>.

Lee, JoonHo. 2026. “The Geometry of Informativeness: A Bregman-Centroid
Diagnostic for Partial Pooling in Multisite Trials.” Unpublished
manuscript.

Rickles, Jordan, Margaret Clements, Iliana Brodziak de los Reyes, Mark
Lachowicz, Shuqiong Lin, and Jessica Heppen. 2024. “A Multisite
Randomized Study of an Online Learning Approach to High School Credit
Recovery: Effects on Student Experiences and Proximal Outcomes.”
*Journal of Research on Educational Effectiveness* 17 (3): 467–90.
<https://doi.org/10.1080/19345747.2023.2198524>.

Rubin, Donald B. 1981. “Estimation in Parallel Randomized Experiments.”
*Journal of Educational Statistics* 6 (4): 377–401.
<https://doi.org/10.3102/10769986006004377>.
