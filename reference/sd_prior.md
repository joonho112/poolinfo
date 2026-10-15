# Prior distribution of the between-site standard deviation

Specifies the prior of `sigma_tau` for
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md):
the flat prior of Lee (2026), which is the default, or a half-normal,
half-t, half-Cauchy or uniform prior.

## Usage

``` r
sd_prior(
  type = c("flat", "half_normal", "half_t", "half_cauchy", "uniform"),
  scale = NULL,
  df = NULL,
  upper = NULL
)
```

## Arguments

- type:

  Family of the prior: one of `"flat"` (the default), `"half_normal"`,
  `"half_t"`, `"half_cauchy"` or `"uniform"`.

- scale:

  Scale of a half-normal, half-t or half-Cauchy prior, a positive number
  in the units of the estimates. It has no default, because a sensible
  scale depends on the scale of the effects.

- df:

  Degrees of freedom of a half-t prior, a positive number.

- upper:

  Upper limit of a uniform prior, a positive number.

## Value

An object of class `sd_prior`, a list with the components

- `type`: the family.

- `scale`, `df` or `upper`: the parameters that the family uses.

- `tail`: the power \\k\\ with which the density falls off in its upper
  tail, which decides which posterior moments exist: 0 for the flat
  prior, `df + 1` for the half-t, 2 for the half-Cauchy and `Inf` for
  the half-normal and uniform priors.

## Details

Proper SD priors are an experimental extension beyond Lee (2026); their
interfaces may change in a later 0.x release.

The flat prior is uniform on the mean of the effects and on `sigma_tau`
(not on `sigma_tau^2`). Under it the posterior is proper only with at
least 3 sites, and the posterior mean of \\\sigma\_\tau^{q}\\ exists
only with more than \\q + 2\\ sites; so with 3 or 4 sites the posterior
mean of `sigma_tau^2` does not exist, although those of the bounded
diagnostics do. With a mean model of \\p\\ columns, the condition is \\J
\> q + p + 1\\.

A prior whose density falls off like \\\sigma^{-k}\\ lowers the bound by
\\k\\: the moment of order \\q\\ exists when \\J \> q + p + 1 - k\\,
with \\k = \nu + 1\\ for a half-t prior with \\\nu\\ degrees of freedom
and \\k = 2\\ for the half-Cauchy. Every moment exists under the
half-normal and the uniform priors. Every prior other than the flat one
needs at least \\p + 1\\ sites.

With few sites the flat prior can matter, as for any prior on a variance
parameter (Gelman, 2006);
[`prior_sensitivity()`](https://joonho112.github.io/poolinfo/reference/prior_sensitivity.md)
compares the summaries under several priors.

## References

Gelman, A. (2006). Prior distributions for variance parameters in
hierarchical models (comment on article by Browne and Draper). *Bayesian
Analysis, 1*(3), 515-534.
[doi:10.1214/06-BA117A](https://doi.org/10.1214/06-BA117A)

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`prior_sensitivity()`](https://joonho112.github.io/poolinfo/reference/prior_sensitivity.md)
to compare the summaries under several priors, and
[`vignette("posterior-summaries")`](https://joonho112.github.io/poolinfo/articles/posterior-summaries.md)
for the choice of prior.

## Examples

``` r
sd_prior()
#> Prior: flat on the mean and on sigma_tau
sd_prior("half_normal", scale = 0.5)
#> Prior: half-normal on sigma_tau with scale 0.5, flat on the mean
sd_prior("half_t", scale = 0.5, df = 3)
#> Prior: half-t on sigma_tau with scale 0.5 and 3 degrees of freedom, flat on the mean
```
