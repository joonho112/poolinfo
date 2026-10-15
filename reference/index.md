# Package index

## Posterior summaries

The posterior distribution of the diagnostics given the site estimates
and their standard errors, and what can be read from it. Proper SD
priors and mean models with covariates are experimental extensions.

- [`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  [`format(`*`<reliability_posterior>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  [`print(`*`<reliability_posterior>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  [`summary(`*`<reliability_posterior>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  [`format(`*`<summary.reliability_posterior>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  [`print(`*`<summary.reliability_posterior>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  [`update(`*`<reliability_posterior>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  [`as.data.frame(`*`<reliability_posterior>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  : Posterior summaries of the reliability profile
- [`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md)
  : Prior distribution of the between-site standard deviation
- [`posterior_control()`](https://joonho112.github.io/poolinfo/reference/posterior_control.md)
  : Settings of the posterior computation
- [`profile_summary()`](https://joonho112.github.io/poolinfo/reference/profile_summary.md)
  : Summaries of a profile as a data frame
- [`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md)
  : Site table of a profile
- [`prob_above()`](https://joonho112.github.io/poolinfo/reference/prob_above.md)
  : Posterior probability that a diagnostic exceeds a value
- [`plugin_gap()`](https://joonho112.github.io/poolinfo/reference/plugin_gap.md)
  : Posterior means against plug-in values
- [`functional_density()`](https://joonho112.github.io/poolinfo/reference/functional_density.md)
  : Posterior density of a diagnostic
- [`profile_draws()`](https://joonho112.github.io/poolinfo/reference/profile_draws.md)
  : Draws from the posterior
- [`check_numerics()`](https://joonho112.github.io/poolinfo/reference/check_numerics.md)
  : Numerical checks of a posterior fit

## Profiles at a fixed between-site SD

The reliability profile when the between-site standard deviation is
given, for example from a REML fit or a planning value.

- [`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  [`update(`*`<reliability_profile>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  [`format(`*`<reliability_profile>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  [`print(`*`<reliability_profile>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  [`summary(`*`<reliability_profile>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  [`format(`*`<summary.reliability_profile>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  [`print(`*`<summary.reliability_profile>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  [`as.data.frame(`*`<reliability_profile>`*`)`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  : Reliability profile at a fixed between-site SD

## Fits from other packages

Diagnostics from posterior draws, from a confidence interval for the
between-site variance, and from fits of metafor, ebnm and multisiteDGP.
External-draw summaries and confidence-interval mapping are experimental
extensions whose interfaces may change in a later 0.x release.

- [`reliability_draws()`](https://joonho112.github.io/poolinfo/reference/reliability_draws.md)
  : Diagnostics from posterior draws of the between-site SD
- [`map_tau2_interval()`](https://joonho112.github.io/poolinfo/reference/map_tau2_interval.md)
  : Confidence interval for the between-site variance mapped through the
  diagnostics
- [`reliability_profile(`*`<rma.uni>`*`)`](https://joonho112.github.io/poolinfo/reference/metafor-methods.md)
  [`reliability_profile(`*`<rma>`*`)`](https://joonho112.github.io/poolinfo/reference/metafor-methods.md)
  : Profiles and posteriors of metafor fits
- [`reliability_profile(`*`<ebnm>`*`)`](https://joonho112.github.io/poolinfo/reference/ebnm-methods.md)
  : Profiles of ebnm fits with a normal prior
- [`reliability_profile(`*`<multisitedgp_data>`*`)`](https://joonho112.github.io/poolinfo/reference/multisitedgp-methods.md)
  : Profiles and posteriors of simulated multisite data

## Plots

Every plot is a ggplot object, or a patchwork object for several panels.

### Analysis

- [`plot_profile()`](https://joonho112.github.io/poolinfo/reference/plot_profile.md)
  : Reliability profile plot
- [`plot_pooling()`](https://joonho112.github.io/poolinfo/reference/plot_pooling.md)
  : Pooling plot
- [`plot_functionals()`](https://joonho112.github.io/poolinfo/reference/plot_functionals.md)
  : Plot of the posterior density of a diagnostic
- [`plot_plane()`](https://joonho112.github.io/poolinfo/reference/plot_plane.md)
  : Informativeness and logit variance plane
- [`plot_dashboard()`](https://joonho112.github.io/poolinfo/reference/plot_dashboard.md)
  : Dashboard of a posterior fit
- [`plot_compare()`](https://joonho112.github.io/poolinfo/reference/plot_compare.md)
  : Comparison plot
- [`plot_sigma_scale()`](https://joonho112.github.io/poolinfo/reference/plot_sigma_scale.md)
  : Posterior of sigma_tau on the scale of the informativeness
- [`plot(`*`<reliability_profile>`*`)`](https://joonho112.github.io/poolinfo/reference/plot.poolinfo.md)
  [`plot(`*`<reliability_posterior>`*`)`](https://joonho112.github.io/poolinfo/reference/plot.poolinfo.md)
  [`plot(`*`<reliability_draws>`*`)`](https://joonho112.github.io/poolinfo/reference/plot.poolinfo.md)
  [`plot(`*`<reliability_design>`*`)`](https://joonho112.github.io/poolinfo/reference/plot.poolinfo.md)
  : Plots of poolinfo objects

### Diagnostics

- [`plot_gap_path()`](https://joonho112.github.io/poolinfo/reference/plot_gap_path.md)
  : Dispersion and centroid gap along the informativeness
- [`plot_concentration()`](https://joonho112.github.io/poolinfo/reference/plot_concentration.md)
  : Concentration of the reliability
- [`plot_contributions()`](https://joonho112.github.io/poolinfo/reference/plot_contributions.md)
  : Sites' terms of the dispersion

### Geometry

- [`plot_geometry()`](https://joonho112.github.io/poolinfo/reference/plot_geometry.md)
  : Geometry of the Kullback-Leibler centroid
- [`plot_decomposition()`](https://joonho112.github.io/poolinfo/reference/plot_decomposition.md)
  : Decomposition of the average divergence

### Design

- [`plot_design()`](https://joonho112.github.io/poolinfo/reference/plot_design.md)
  : Design plots

## Comparing analyses

Several outcomes, scales, subsets or priors side by side.

- [`compare_profiles()`](https://joonho112.github.io/poolinfo/reference/compare_profiles.md)
  [`plot(`*`<profile_comparison>`*`)`](https://joonho112.github.io/poolinfo/reference/compare_profiles.md)
  : Diagnostics of several analyses side by side
- [`prior_sensitivity()`](https://joonho112.github.io/poolinfo/reference/prior_sensitivity.md)
  : Posterior summaries under several priors

## Reporting

Tables and sentences for a report, with the conditioning they rest on.

- [`report_table()`](https://joonho112.github.io/poolinfo/reference/report_table.md)
  : Table of the diagnostics for a report
- [`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md)
  : Sentences for a report

## Planning

Sampling variances, profiles and allocations of a trial before any data
are collected.

- [`neyman_variance()`](https://joonho112.github.io/poolinfo/reference/neyman_variance.md)
  : Sampling variances of a planned trial
- [`reliability_design()`](https://joonho112.github.io/poolinfo/reference/reliability_design.md)
  : Reliability profile of a planned allocation
- [`balanced_allocation()`](https://joonho112.github.io/poolinfo/reference/balanced_allocation.md)
  : Site allocations under size and budget constraints
- [`solve_design()`](https://joonho112.github.io/poolinfo/reference/solve_design.md)
  : Sizes that reach a target informativeness
- [`imbalance_ratio()`](https://joonho112.github.io/poolinfo/reference/imbalance_ratio.md)
  : Ratio of the geometric to the arithmetic mean site size

## Diagnostics of a profile

The informativeness, the mean reliability, the dispersion and the
concentration count of a set of site reliabilities, and the quantities
related to them.

- [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md)
  : Site reliabilities
- [`kl_centroid()`](https://joonho112.github.io/poolinfo/reference/kl_centroid.md)
  : Informativeness of a reliability profile
- [`reliability_centroids()`](https://joonho112.github.io/poolinfo/reference/reliability_centroids.md)
  : Three centroids of a reliability profile
- [`kl_dispersion()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md)
  [`kl_contributions()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md)
  : Dispersion of a reliability profile
- [`concentration_count()`](https://joonho112.github.io/poolinfo/reference/concentration_count.md)
  : Concentration count of a reliability profile
- [`logit_moments()`](https://joonho112.github.io/poolinfo/reference/logit_moments.md)
  : Logit variance and third absolute moment
- [`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md)
  : Second-order terms of the dispersion and of the centroid gap
- [`typical_variance()`](https://joonho112.github.io/poolinfo/reference/typical_variance.md)
  : Typical sampling variance of Higgins and Thompson
- [`benchmark_risk()`](https://joonho112.github.io/poolinfo/reference/benchmark_risk.md)
  : Benchmark risk of a reliability profile
- [`bregman_centroid()`](https://joonho112.github.io/poolinfo/reference/bregman_centroid.md)
  [`bregman_dispersion()`](https://joonho112.github.io/poolinfo/reference/bregman_centroid.md)
  : Bregman centroids and dispersions of a reliability profile
- [`profile_path()`](https://joonho112.github.io/poolinfo/reference/profile_path.md)
  : Dispersion and centroid gap along a path of informativeness
- [`shape_profile()`](https://joonho112.github.io/poolinfo/reference/shape_profile.md)
  : Sampling variances with a given informativeness and logit variance

## Data sets

- [`ocrs`](https://joonho112.github.io/poolinfo/reference/ocrs.md) :
  Online Credit Recovery Study
- [`eight_schools`](https://joonho112.github.io/poolinfo/reference/eight_schools.md)
  : Eight schools
- [`example_profiles`](https://joonho112.github.io/poolinfo/reference/example_profiles.md)
  : Example reliability profiles

## Package

- [`poolinfo`](https://joonho112.github.io/poolinfo/reference/poolinfo-package.md)
  [`poolinfo-package`](https://joonho112.github.io/poolinfo/reference/poolinfo-package.md)
  : poolinfo: Reliability-Profile Diagnostics for Partial Pooling
