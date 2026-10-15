# Settings of the posterior computation

Chooses the engine that integrates over the between-site standard
deviation and sets its tolerances.

## Usage

``` r
posterior_control(
  engine = c("nodes", "adaptive"),
  rel_tol = 1e-10,
  abs_tol = 1e-12,
  chunk_sites = 2000
)
```

## Arguments

- engine:

  The engine: `"nodes"` (the default) builds one quadrature rule and
  reuses it for every summary of a fit; `"adaptive"` computes each
  posterior mean, probability and quantile with
  [`stats::integrate()`](https://rdrr.io/r/stats/integrate.html), as Lee
  (2026, Appendix D.3) describes.

- rel_tol, abs_tol:

  Relative and absolute tolerances of
  [`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) for the
  adaptive engine, positive numbers. The absolute tolerance applies to
  integrands scaled to the size of the posterior and of the integrated
  quantity, so it does not depend on the units of the estimates.

- chunk_sites:

  Number of sites computed at a time for the site summaries, a positive
  integer. Smaller values use less memory.

## Value

A list of class `posterior_control` with the components

- `engine`: `"nodes"` or `"adaptive"`.

- `rel_tol`, `abs_tol`: the tolerances of the adaptive engine.

- `chunk_sites`: the number of sites computed at a time.

## Details

Both engines integrate the same posterior, and their summaries agree to
many more digits than are printed. The node rule is the faster: one rule
of Gauss-Legendre panels serves every summary of a fit, while the
adaptive engine integrates each summary, and each site, separately. The
adaptive engine takes the posterior standard deviations, and the
intervals of the site effects and of the coefficients, from a node rule.
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
describes the rule, its checks and the comparison of the engines. With a
prior whose tail is very heavy, such as a half-t prior with far fewer
than one degree of freedom, a moment can exist and still have a tail too
slow for the rule to bound. Both engines then stop with an error that
names the moment. A moment that does not exist is `NA` instead.
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
gives the limits.

## See also

[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md),
whose argument `control` takes these settings, and
[`check_numerics()`](https://joonho112.github.io/poolinfo/reference/check_numerics.md)
for the checks of a fit.

## Examples

``` r
posterior_control()
#> Posterior computation: nodes engine
posterior_control(engine = "adaptive")
#> Posterior computation: adaptive engine, rel_tol 1e-10, abs_tol 1e-12
```
