# Plots of poolinfo objects

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows one of
the package's plots of a profile, a posterior fit, draws or a design.
`type` chooses the plot, and most types call the plot function of the
same name, so `plot(fit, type = "pooling")` is `plot_pooling(fit)`.

## Usage

``` r
# S3 method for class 'reliability_profile'
plot(
  x,
  type = c("profile", "logit", "geometry", "decomposition", "gap", "concentration",
    "contributions"),
  ...
)

# S3 method for class 'reliability_posterior'
plot(
  x,
  type = c("dashboard", "profile", "logit", "pooling", "functionals", "plane", "sigma",
    "gap", "concentration", "contributions"),
  ...
)

# S3 method for class 'reliability_draws'
plot(
  x,
  type = c("dashboard", "profile", "logit", "pooling", "functionals", "plane", "sigma",
    "gap", "concentration", "contributions"),
  ...
)

# S3 method for class 'reliability_design'
plot(x, type = c("imbalance", "sites", "target", "profile"), ...)
```

## Arguments

- x:

  A `reliability_profile`, `reliability_posterior`, `reliability_draws`
  or `reliability_design` object.

- type:

  The plot:

  - for a profile, `"profile"` (the default), `"logit"`, `"geometry"`,
    `"decomposition"`, `"gap"`, `"concentration"` or `"contributions"`.

  - for a posterior fit or draws, `"dashboard"` (the default),
    `"profile"`, `"logit"`, `"pooling"`, `"functionals"`, `"plane"`,
    `"sigma"`, `"gap"`, `"concentration"` or `"contributions"`.

  - for a design, `"imbalance"` (the default), `"sites"`, `"target"` or
    `"profile"`.

  `"logit"` is
  [`plot_profile()`](https://joonho112.github.io/poolinfo/reference/plot_profile.md)
  with `type = "logit"`, `"sigma"` is
  [`plot_sigma_scale()`](https://joonho112.github.io/poolinfo/reference/plot_sigma_scale.md),
  `"gap"` is
  [`plot_gap_path()`](https://joonho112.github.io/poolinfo/reference/plot_gap_path.md),
  and `"imbalance"`, `"sites"` and `"target"` are the types of
  [`plot_design()`](https://joonho112.github.io/poolinfo/reference/plot_design.md);
  `"profile"` of a design shows the planned profile.

- ...:

  Further arguments of the plot function.

## Value

A ggplot or patchwork object.

## See also

[`plot_profile()`](https://joonho112.github.io/poolinfo/reference/plot_profile.md),
[`plot_dashboard()`](https://joonho112.github.io/poolinfo/reference/plot_dashboard.md),
[`plot_design()`](https://joonho112.github.io/poolinfo/reference/plot_design.md),
and
[`compare_profiles()`](https://joonho112.github.io/poolinfo/reference/compare_profiles.md)
for the plots of a comparison.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
plot(fit)

plot(fit, type = "functionals")
```
