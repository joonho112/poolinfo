# poolinfo

A multisite trial, a meta-analysis or a small area estimation problem
gives an estimate and a standard error for each site. A hierarchical
model shrinks each estimate toward the overall mean, and the weight that
a site’s posterior mean keeps on its own estimate is the site’s
reliability. poolinfo summarizes the reliabilities of all the sites at
once: their informativeness (the reliability at the geometric mean of
the sampling variances), their mean, their dispersion, and a count of
how the total reliability is shared among the sites. It computes the
posterior distribution of each from the estimates and their standard
errors, and it plots them. The summaries are worth checking before you
rank sites, report site effects or study how an effect varies across
sites, because they measure how much of each site’s own estimate
survives the pooling. The diagnostics are those of Lee (2026).

## Installation

You can install poolinfo from GitHub with:

``` r

# install.packages("remotes")
remotes::install_github("joonho112/poolinfo", build_vignettes = TRUE)
```

The vignettes are also on the [package
website](https://joonho112.github.io/poolinfo/).

## Example

The data set `ocrs` has the estimated effects of an online credit
recovery course on course passing in the 32 school-by-subject sites of
the Online Credit Recovery Study ([Rickles et al.,
2024](https://doi.org/10.1080/19345747.2023.2198524)), with their
standard errors.
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
computes the posterior distribution of the diagnostics by numerical
integration over the between-site standard deviation, with a flat prior
on it and on the mean.

``` r

library(poolinfo)
fit <- reliability_posterior(estimate, se = se, site = site, data = ocrs,
                             scale_label = "risk difference")
fit
#> Reliability profile of 32 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Scale of the effects: risk difference
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

The informativeness is 0.67, with a 95% posterior interval from 0.50 to
0.81, so a site whose sampling variance is the geometric mean of the 32
keeps about 67% of its own estimate in its posterior mean. The mean
reliability, 0.66, is close to it, the dispersion is H = 0.025, and the
concentration count is 31.2 of 32, so the total reliability is shared
almost evenly among the sites. Every summary is conditional on the
standard errors, which the model treats as known.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows the site
estimates and their posterior means, the sorted reliabilities, the
posterior density of the informativeness and the plane of the
informativeness and the logit variance.

``` r

plot(fit, scale = 100)
```

![Four panels: the site estimates and their posterior means with
intervals; the sorted site reliabilities with the informativeness and
its band; the posterior density of the informativeness; and the
informativeness and logit variance
plane.](reference/figures/README-dashboard-1.png)

## What else is in the package

- Diagnostics of a profile at a fixed between-site standard deviation:
  [`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md),
  [`kl_centroid()`](https://joonho112.github.io/poolinfo/reference/kl_centroid.md),
  [`kl_dispersion()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md),
  [`concentration_count()`](https://joonho112.github.io/poolinfo/reference/concentration_count.md)
  and related functions.
- Diagnostics from fits of other packages:
  [`reliability_draws()`](https://joonho112.github.io/poolinfo/reference/reliability_draws.md)
  for posterior draws (from brms, for example), and methods for metafor,
  ebnm and multisiteDGP objects;
  [`map_tau2_interval()`](https://joonho112.github.io/poolinfo/reference/map_tau2_interval.md)
  maps a confidence interval for the between-site variance.
- Planning:
  [`reliability_design()`](https://joonho112.github.io/poolinfo/reference/reliability_design.md),
  [`balanced_allocation()`](https://joonho112.github.io/poolinfo/reference/balanced_allocation.md)
  and
  [`solve_design()`](https://joonho112.github.io/poolinfo/reference/solve_design.md)
  give the profile of a planned trial and allocations of a total or a
  budget to the sites.
- Comparing and reporting:
  [`compare_profiles()`](https://joonho112.github.io/poolinfo/reference/compare_profiles.md),
  [`prior_sensitivity()`](https://joonho112.github.io/poolinfo/reference/prior_sensitivity.md),
  [`report_table()`](https://joonho112.github.io/poolinfo/reference/report_table.md)
  and
  [`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md).

[`vignette("poolinfo")`](https://joonho112.github.io/poolinfo/articles/poolinfo.md)
walks through the example above in full.

## Citation

To cite poolinfo, use `citation("poolinfo")`. The diagnostics are those
of

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Funding

This research was supported by the Institute of Education Sciences, U.S.
Department of Education, through Grant R305D240078 to the University of
Alabama. The opinions expressed are those of the author and do not
represent views of the Institute or the U.S. Department of Education.
