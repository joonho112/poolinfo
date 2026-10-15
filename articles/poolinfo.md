# Getting started with poolinfo

In a multisite trial each site gives an estimate of the intervention’s
effect, with a standard error. A hierarchical model shrinks each
estimate toward the overall mean, a noisy site’s estimate more than a
precise one’s. poolinfo summarizes this shrinkage across all the sites
of a trial. The summaries are worth a look before you rank the sites,
report their effects or study how the effect varies among them, because
they say how much of each site’s own estimate survives in the model’s
results, and they bear on which estimator suits such a goal ([Lee
2026](#ref-lee2026)). In this vignette we analyze one outcome of a real
trial, from the site estimates to the report.

``` r

library(poolinfo)
```

## The data

The Online Credit Recovery Study randomized students who had failed
Algebra 1 or English 9 in their first year of high school to one of two
summer courses at their school, an online course or the school’s usual
teacher-directed course ([Rickles et al. 2024](#ref-rickles2024)). The
data set `ocrs` has one row for each of the 32 school-by-subject sites
in which students were assigned to both courses.

``` r

head(ocrs[, c("site", "subject", "n", "estimate", "se")])
#>   site   subject  n    estimate         se
#> 1  S01 Algebra 1 57  0.20526372 0.13001694
#> 2  S02 Algebra 1 47 -0.27638587 0.10163046
#> 3  S03 Algebra 1 58  0.03845466 0.09179922
#> 4  S04 Algebra 1 39 -0.12344656 0.16209038
#> 5  S05 Algebra 1 41 -0.01560407 0.15503087
#> 6  S06 Algebra 1 38  0.06041924 0.17049175
```

Each `estimate` is the effect of assignment to the online course on
passing the course, as a difference in proportions, and `se` is its
standard error. A negative estimate means that, in that site, a smaller
proportion of the students assigned to the online course passed. Within
a site, the estimate comes from a regression of passing on the
assignment indicator and the randomization blocks, so it weights the
blocks differently from a comparison of all the site’s students.
[`?ocrs`](https://joonho112.github.io/poolinfo/reference/ocrs.md) gives
the details, and
[`vignette("multisite-trial")`](https://joonho112.github.io/poolinfo/articles/multisite-trial.md)
compares this contrast with other first-stage analyses.

## Fitting the posterior

[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
takes the estimates and their standard errors. We also record the scale
of the effects and the estimand at each site, which the printed output
and the report sentences repeat.

``` r

fit <- reliability_posterior(
  estimate, se = se, site = site, data = ocrs,
  scale_label = "risk difference",
  target_label = "the block-adjusted intention-to-treat contrast for course passing"
)
fit
#> Reliability profile of 32 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Scale of the effects: risk difference
#> Site estimates: the block-adjusted intention-to-treat contrast for course
#>                 passing
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#>                              mean     95% interval
#>   informativeness I          0.671    0.497    0.812
#>   mean reliability           0.663    0.498    0.801
#>   dispersion H               0.0253   0.0185   0.0289
#>   concentration count J_eff  31.2     30.3     31.7     of 32
#>   logit variance V           0.237    (fixed by the standard errors)
#>   sigma_tau                  0.188    0.127    0.266
#>   effect mean                -0.130   -0.211   -0.0498
#> 
#>   Pr(I > 0.5)                0.973
#>   I at E(sigma_tau^2)        0.691
```

We use flat priors on the mean and on the between-site standard
deviation $`\sigma_\tau`$, as Lee ([2026](#ref-lee2026)) does, and the
posterior comes from numerical integration over $`\sigma_\tau`$ in a
fraction of a second. See
[`vignette("posterior-summaries")`](https://joonho112.github.io/poolinfo/articles/posterior-summaries.md)
for other priors and how much they matter. The last two lines of the
header state an assumption and a definition. The model treats the
standard errors as exact, although they are estimates too, so the
posterior intervals leave out the uncertainty in the standard errors
themselves, and a different first-stage analysis of the same trial could
give different diagnostics. And each site counts equally in the
diagnostics, whatever its size, so a site of 24 students weighs as much
as one of 104.

## Site reliabilities

Write $`\hat\tau_j`$ for site $`j`$’s estimate and $`v_j`$ for its
sampling variance, the squared standard error. Given $`\sigma_\tau`$ and
the mean $`\mu`$ of the effects, the posterior mean of site $`j`$’s
effect is

``` math
S_j \hat\tau_j + (1 - S_j)\mu, \qquad
S_j = \frac{\sigma_\tau^2}{\sigma_\tau^2 + v_j}.
```

$`S_j`$ is the weight that the normal model puts on the site’s own
estimate, and it is also the reliability of the site’s estimate, the
proportion of the estimate’s variance that comes from the true effects.
We call the reliabilities $`S_1, \dots, S_J`$ together the trial’s
reliability profile. Lee ([2026](#ref-lee2026)) writes $`\sigma^2`$ and
$`\tau`$ for what this package calls $`\sigma_\tau^2`$ and $`\mu`$.

## The informativeness

The first summary of the profile is the *informativeness*,

``` math
I = \frac{\sigma_\tau^2}{\sigma_\tau^2 + \mathrm{GM}(v)},
```

the reliability of a site whose sampling variance is the geometric mean
$`\mathrm{GM}(v)`$ of the $`v_j`$. It is also the mean of the
reliabilities on the logit scale ([Lee 2026](#ref-lee2026)), and Lee et
al. ([2025](#ref-lee2025)) use it to organize comparisons of prior
models and posterior summaries for site effects.

For course passing, the posterior mean of $`I`$ is 0.67, with a 95%
posterior interval from 0.50 to 0.81. The geometric mean of the 32
standard errors is 0.128, and the posterior mean of $`\sigma_\tau`$ is
0.19, so $`\sigma_\tau^2`$ is about 2.2 times $`\mathrm{GM}(v)`$. A site
with the geometric-mean standard error keeps about 67% of its own
estimate, which means that the model moves its posterior mean about 33%
of the way toward the mean effect. The posterior probability that $`I`$
exceeds 1/2 is 0.97. At $`I = 1/2`$, $`\sigma_\tau^2`$ equals
$`\mathrm{GM}(v)`$ and a site with the geometric-mean standard error is
pulled halfway.

Meta-analysts will recognize the form of the $`I^2`$ statistic ([Higgins
and Thompson 2002](#ref-higgins2002)), which puts the “typical” sampling
variance where $`I`$ puts the geometric mean. For these estimates,
metafor’s REML fit gives an $`I^2`$ of 69%, and its estimate of
$`\sigma_\tau`$, 0.179, is a little below the posterior mean.

## How unequal the reliabilities are

The other three summaries describe how the reliabilities differ from
site to site, and all three depend on the logit variance $`V`$, the
variance of the logits of the reliabilities (with divisor $`J`$).
Because $`\mathrm{logit}\, S_j = \log\sigma_\tau^2 - \log v_j`$, $`V`$
is also the variance of $`\log v_j`$, so the standard errors fix it.
Here $`V = 0.24`$, so the standard errors are typically within a factor
of 1.28 of their geometric mean, although they range from 0.077 to
0.218.

The *mean reliability* $`\bar S`$ is the ordinary average of the
$`S_j`$. With $`\sigma_\tau`$ and $`\mu`$ known, the average posterior
variance of the site effects is $`\sigma_\tau^2(1 - \bar S)`$. In a
trial whose sites all had the geometric-mean standard error, every
$`S_j`$ would equal $`I`$, so $`\bar S > I`$ exactly when the site
effects are estimated more precisely, on average, than in that trial of
equal precision. For passing, $`\bar S`$ is 0.66, a little below $`I`$.

The *dispersion* $`H`$ measures how far the reliabilities are from one
another. Treat each $`S_j`$ as the success probability of a Bernoulli
distribution, a device for measuring how far apart reliabilities are
rather than a model for the data. Among constant values $`p`$, the
average Kullback-Leibler divergence
$`J^{-1}\sum_j \mathrm{KL}\{\mathrm{Ber}(p) \,\|\, \mathrm{Ber}(S_j)\}`$
is smallest at $`p = I`$, and $`H`$ is that smallest value. It is 0 only
when all the reliabilities are equal. When they are close to one
another, the leading terms in $`V`$ are

``` math
H \approx \tfrac12 I(1 - I)V, \qquad
\bar S - I \approx \tfrac12 I(1 - I)(1 - 2I)V.
```

The mean reliability is therefore usually below $`I`$ when $`I > 1/2`$
and above it when $`I < 1/2`$, but a skewed profile, or a value of $`I`$
near 1/2, can reverse the sign.
[`vignette("reliability-profiles")`](https://joonho112.github.io/poolinfo/articles/reliability-profiles.md)
shows when. For passing, the leading term of $`H`$ at the posterior mean
of $`I`$ is 0.026, close to $`H`$’s posterior mean of 0.025.

The *concentration count*

``` math
J_{\mathrm{eff}} = \frac{(\sum_j S_j)^2}{\sum_j S_j^2}
```

describes how the total reliability $`\sum_j S_j`$ is shared among the
sites. It equals $`J`$ when every site has the same reliability,
whatever that reliability is, and it approaches 1 as one site comes to
carry nearly all of it. Its formula is Kish’s, with the reliabilities as
weights, but it is not a sample size. It measures only how evenly the
total is shared, so sites that all have a reliability of 0.1 give $`J`$,
and so do sites that all have 0.9. For passing it is 31.2 of 32, so the
total is shared almost evenly.

## The posterior in more detail

[`summary()`](https://rdrr.io/r/base/summary.html) adds the posterior
standard deviations and medians, the posterior of $`\sigma_\tau`$ and of
the mean, and each diagnostic’s value at the posterior mean of
$`\sigma_\tau^2`$.

``` r

summary(fit)
#> Reliability profile of 32 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Scale of the effects: risk difference
#> Site estimates: the block-adjusted intention-to-treat contrast for course
#>                 passing
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#> Diagnostics
#>                              mean      sd        2.5%      50%       97.5%
#>   informativeness I          0.671     0.081     0.497     0.677     0.812
#>   mean reliability           0.663     0.078     0.498     0.668     0.801
#>   harmonic centroid          0.644     0.084     0.467     0.650     0.793
#>   dispersion H               0.0253    0.0029    0.0185    0.0258    0.0289
#>   concentration count J_eff  31.2      0.4       30.3      31.2      31.7
#>   logit variance V           0.237     (fixed by the standard errors)
#> 
#> Between-site SD and mean
#>                              mean      sd        2.5%      50%       97.5%
#>   sigma_tau                  0.188     0.0354    0.127     0.185     0.266
#>   sigma_tau^2                0.0365    0.0141    0.0162    0.0342    0.0706
#>   effect mean                -0.130    0.0409    -0.211    -0.131    -0.0498
#> 
#> Pr(I > 0.5) = 0.973
#> 
#> Posterior means and values at the posterior mean of sigma_tau^2
#>                                      posterior at E(sigma_tau^2)
#>   informativeness I                  0.671     0.691
#>   mean reliability                   0.663     0.682
#>   dispersion H                       0.0253    0.0253
#>   concentration count J_eff          31.2      31.3
#> 
#> H is largest, at 0.0289, where sigma_tau = 0.130 (I = 0.509).
#> At that sigma_tau the mean reliability equals I.
#> Integration over sigma_tau: a rule of 800 nodes.
#> All numerical checks are within their limits.
```

The first block also gives the harmonic centroid, the harmonic mean of
the $`S_j`$, which equals $`\sigma_\tau^2/(\sigma_\tau^2 + \bar v)`$
with $`\bar v`$ the arithmetic mean of the $`v_j`$, so it is never above
$`I`$ or $`\bar S`$. To use the numbers in your own code,
[`profile_summary()`](https://joonho112.github.io/poolinfo/reference/profile_summary.md)
returns them as a data frame, one row per quantity.

``` r

print(profile_summary(fit), digits = 3)
#>      quantity    mean      sd   lower   upper
#> 1           I  0.6709 0.08079  0.4974  0.8122
#> 2        mean  0.6635 0.07775  0.4980  0.8013
#> 3    harmonic  0.6445 0.08362  0.4672  0.7930
#> 4           H  0.0253 0.00286  0.0185  0.0289
#> 5       J_eff 31.1704 0.35109 30.3456 31.7004
#> 6           V  0.2369 0.00000  0.2369  0.2369
#> 7   sigma_tau  0.1879 0.03538  0.1272  0.2658
#> 8  sigma_tau2  0.0365 0.01411  0.0162  0.0706
#> 9 effect_mean -0.1305 0.04089 -0.2111 -0.0498
```

The mean effect has posterior mean -13.0 percentage points, with a 95%
posterior interval from -21.1 to -5.0. It is the mean of the
distribution of site effects, a mean over sites rather than students, so
it is not the effect for the average student. Because its interval lies
below 0, we conclude that assignment to the online course lowered
passing in the average site.

In the last block we compare each posterior mean with the diagnostic’s
value at the posterior mean of $`\sigma_\tau^2`$, the number that a
plug-in of that estimate would report. $`I`$ and every $`S_j`$ are
concave functions of $`\sigma_\tau^2`$, so $`I`$ and $`\bar S`$ at the
posterior mean of $`\sigma_\tau^2`$ are never below their posterior
means ([Lee 2026](#ref-lee2026), Proposition 3). Here $`I`$ is 0.691
against 0.671. $`H`$ and $`J_{\mathrm{eff}}`$ have no such order. A
plug-in of another estimate of $`\sigma_\tau`$ can fall on either side
of the posterior mean: at metafor’s REML estimate, $`I`$ is 0.663.

## Plots

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows four
panels: the estimates beside the posterior means of the site effects,
the sorted reliabilities, the posterior density of $`I`$, and the
analysis’s position in the plane of $`I`$ and $`V`$. We multiply the
effects by 100 to show them in percentage points.

``` r

plot(fit, scale = 100)
```

![Four panels: the site estimates and their posterior means with
intervals; the sorted site reliabilities with the informativeness and
its band; the posterior density of the informativeness; and the
informativeness and logit variance
plane.](poolinfo_files/figure-html/dashboard-1.png)

In panel (a) we see the posterior means pulled toward the mean effect,
most visibly for the sites with wide intervals. Panel (b) shows the
reliability profile, from 0.42 to 0.84. In panel (d) the point lies
right of 1/2, where the leading term of $`\bar S - I`$ is negative, and
between the contours at 0.02 and 0.03 of $`H`$’s leading term. The plot
is a patchwork of ggplot objects. To change one panel, call its own
function
([`plot_pooling()`](https://joonho112.github.io/poolinfo/reference/plot_pooling.md),
[`plot_profile()`](https://joonho112.github.io/poolinfo/reference/plot_profile.md),
[`plot_functionals()`](https://joonho112.github.io/poolinfo/reference/plot_functionals.md)
or
[`plot_plane()`](https://joonho112.github.io/poolinfo/reference/plot_plane.md)),
which returns a ggplot object with a fuller caption.

## The sites

[`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md)
gives the numbers behind panels (a) and (b), one row per site. We show
the three sites with the largest standard errors and the three with the
smallest, in percentage points:

``` r

sites <- site_table(fit)
sites <- sites[order(sites$se, decreasing = TRUE),
               c("site", "estimate", "se", "reliability", "effect")]
shown <- rbind(head(sites, 3), tail(sites, 3))
points <- c("estimate", "se", "effect")
shown[points] <- round(100 * shown[points], 1)
shown$reliability <- round(shown$reliability, 2)
shown
#>    site estimate   se reliability effect
#> 24  S24    -41.0 21.8        0.42  -24.8
#> 14  S14    -22.2 20.1        0.46  -17.3
#> 18  S18    -23.1 18.8        0.49  -18.0
#> 22  S22    -33.1  9.3        0.79  -28.9
#> 3   S03      3.8  9.2        0.80    0.4
#> 19  S19    -45.3  7.7        0.84  -40.3
```

Each posterior mean moves toward the mean effect by about $`1 - S_j`$ of
the distance between the site’s estimate and the mean. A noisy site
moves a larger part of the way, but it moves far only if its estimate is
far from the mean. Site S24, with the largest standard error, moves 58%
of the way, from -41.0 to -24.8 points. Site S14 is nearly as noisy and
moves 54% of the way, yet its estimate is close to the mean, so its
posterior mean changes by only 4.9 points.

## Your own data

The same call works for any estimates with standard errors. You can give

- plain vectors, `reliability_posterior(y, se = s)`, or sampling
  variances with `se2`, such as the `vi` column that
  [`metafor::escalc()`](https://wviechtb.github.io/metafor/reference/escalc.html)
  computes;
- a random-effects fit of metafor, `reliability_posterior(rma_fit)`,
  whose moderators become the mean model; `reliability_profile(rma_fit)`
  gives the profile at its REML estimate, and
  [`vignette("meta-analysis")`](https://joonho112.github.io/poolinfo/articles/meta-analysis.md)
  works through an example;
- site covariates for the mean, `mean_model = ~ x`, so that each site is
  shrunk toward its predicted value;
- posterior draws of $`\sigma_\tau`$ from other software, to
  [`reliability_draws()`](https://joonho112.github.io/poolinfo/reference/reliability_draws.md),
  which also reads a brms fit of the form `y | se(sei) ~ 1 + (1 | site)`
  directly.

The model treats the sites as independent. In this trial 8 schools have
a site in each subject, and some students appear in both. We check how
much that matters in
[`vignette("multisite-trial")`](https://joonho112.github.io/poolinfo/articles/multisite-trial.md).

The model also gives every site the same mean and variance of the
effects, whatever its precision. When the smaller sites tend to have
larger effects, say, that common mean and variance can be inadequate
([Lee 2026, sec. 9](#ref-lee2026)). A plot of the estimates against
their standard errors is a first check, although both are noisy:

``` r

ggplot2::ggplot(ocrs, ggplot2::aes(x = se, y = estimate)) +
  ggplot2::geom_point() +
  ggplot2::labs(x = "Standard error", y = "Site estimate") +
  ggplot2::theme_minimal()
```

![The 32 site estimates for course passing plotted against their
standard errors; the points show no
trend.](poolinfo_files/figure-html/estimate-precision-1.png)

In this trial the estimates show no trend with their standard errors (a
correlation of -0.06). A trend would call for a covariate in the mean
model that carries it, or for a model in which the variance of the
effects depends on precision, which poolinfo does not fit.

## Reporting

Lee ([2026](#ref-lee2026)) recommends that a report

1.  state the estimand at each site, its scale, the estimator and the
    standard errors;
2.  give $`I`$ and $`\bar S`$ together, with $`H`$ or $`V`$;
3.  add $`J_{\mathrm{eff}}`$ when the concentration of the reliability
    matters.

[`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md)
writes the numbers as sentences, with the labels and the conditioning on
the standard errors. A report still needs to define the four quantities
and cite the paper.

``` r

writeLines(strwrap(report_text(fit, outcome = "course passing")))
#> The effects are on the risk difference scale.
#> The estimand at each site is the block-adjusted intention-to-treat
#> contrast for course passing.
#> Under the normal working model with flat priors on the mean and the
#> between-site standard deviation, the informativeness of the 32 sites
#> for course passing was 0.67 (95% posterior interval 0.50 to 0.81), and
#> the posterior probability that it exceeded 0.5 was 0.97.
#> At an informativeness of 0.5, the between-site variance equals the
#> geometric mean of the sampling variances.
#> The mean reliability was 0.66 (0.50 to 0.80).
#> The dispersion of the site reliabilities was H = 0.025 (0.019 to
#> 0.029), with logit variance V = 0.24, and the concentration count was
#> 31.2 (at most 32, the number of sites).
#> These summaries are conditional on the supplied standard errors, and
#> each site counts equally in them, whatever its size.
```

[`report_table()`](https://joonho112.github.io/poolinfo/reference/report_table.md)
gives the main numbers as a table, either as a data frame or as markdown
or LaTeX text. In an R Markdown document, a chunk with the option
`results = "asis"` turns the markdown into a table:

``` r

cat(report_table(fit, format = "markdown"))
```

| J | I \[95% interval\] | Mean reliability | H \[95% interval\] | J_eff | Pr(I \> 0.5) |
|:---|:---|:---|:---|:---|:---|
| 32 | .67 \[.50, .81\] | .66 | .025 \[.019, .029\] | 31.2 | .97 |

*Note.* Each estimate is modeled as normal around its site’s effect with
a known sampling variance, and the site effects as independent and
normal across sites. Flat priors on the mean and the between-site
standard deviation. The effects are on the risk difference scale. The
estimand at each site is the block-adjusted intention-to-treat contrast
for course passing. Summaries are conditional on the supplied standard
errors, and each site counts equally in the diagnostics, whatever its
size.

## Other vignettes

- [`vignette("multisite-trial")`](https://joonho112.github.io/poolinfo/articles/multisite-trial.md)
  adds the study’s second outcome, credits earned by the fourth year,
  whose diagnostics differ sharply from those of passing, and shows how
  much the scale of the effects matters.
- [`vignette("reliability-profiles")`](https://joonho112.github.io/poolinfo/articles/reliability-profiles.md)
  works through the geometry behind $`I`$, $`H`$ and
  $`J_{\mathrm{eff}}`$.
- [`vignette("posterior-summaries")`](https://joonho112.github.io/poolinfo/articles/posterior-summaries.md)
  covers the choice of prior, posterior draws from other software and
  the numerical checks.
- [`vignette("meta-analysis")`](https://joonho112.github.io/poolinfo/articles/meta-analysis.md)
  and
  [`vignette("small-area-estimation")`](https://joonho112.github.io/poolinfo/articles/small-area-estimation.md)
  apply the diagnostics to a meta-analysis and to an area-level model.
- [`vignette("planning")`](https://joonho112.github.io/poolinfo/articles/planning.md)
  uses them to plan a trial.
- [`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
  documents the computation.

## References

Higgins, Julian P. T., and Simon G. Thompson. 2002. “Quantifying
Heterogeneity in a Meta-Analysis.” *Statistics in Medicine* 21 (11):
1539–58. <https://doi.org/10.1002/sim.1186>.

Lee, JoonHo. 2026. “The Geometry of Informativeness: A Bregman-Centroid
Diagnostic for Partial Pooling in Multisite Trials.” Unpublished
manuscript.

Lee, JoonHo, Jonathan Che, Sophia Rabe-Hesketh, Avi Feller, and Luke
Miratrix. 2025. “Improving the Estimation of Site-Specific Effects and
Their Distribution in Multisite Trials.” *Journal of Educational and
Behavioral Statistics* 50 (5): 731–64.
<https://doi.org/10.3102/10769986241254286>.

Rickles, Jordan, Margaret Clements, Iliana Brodziak de los Reyes, Mark
Lachowicz, Shuqiong Lin, and Jessica Heppen. 2024. “A Multisite
Randomized Study of an Online Learning Approach to High School Credit
Recovery: Effects on Student Experiences and Proximal Outcomes.”
*Journal of Research on Educational Effectiveness* 17 (3): 467–90.
<https://doi.org/10.1080/19345747.2023.2198524>.
