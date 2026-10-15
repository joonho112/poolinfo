# Example reliability profiles

Five sets of sampling variances for 40 sites whose reliability profiles
at `sigma_tau = 0.2` all have the informativeness 0.3: equal variances,
and four shapes with the logit variance 1.5. The homogeneous, even and
concentrated_precise profiles are those of Figure 2 of Lee (2026).

## Usage

``` r
example_profiles
```

## Format

A data frame with 200 rows and 6 columns:

- profile:

  The profile, a factor: `homogeneous`, `even`, `bimodal`,
  `concentrated_precise` or `concentrated_noisy`; see
  [`shape_profile()`](https://joonho112.github.io/poolinfo/reference/shape_profile.md).

- site:

  Number of the site, 1 to 40.

- se2, se:

  Sampling variance and standard error.

- sigma_tau:

  The between-site standard deviation, 0.2.

- reliability:

  The site's reliability at that `sigma_tau`.

## Source

Built by
[`shape_profile()`](https://joonho112.github.io/poolinfo/reference/shape_profile.md),
with the script `data-raw/example_profiles.R` of the source repository.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
by_profile <- split(example_profiles, example_profiles$profile)
sapply(by_profile, function(d) kl_dispersion(d$reliability))
#>          homogeneous                 even              bimodal 
#>            0.0000000            0.1485641            0.1524161 
#> concentrated_precise   concentrated_noisy 
#>            0.1675842            0.1166457 
```
