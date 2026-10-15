# Sampling variances with a given informativeness and logit variance

Builds sampling variances whose reliability profile at a given
between-site standard deviation has exactly the informativeness `I` and
the logit variance `V`, in one of four shapes. These are the profiles of
the simulation study of Lee (2026, Appendix F).

## Usage

``` r
shape_profile(
  J,
  I,
  V,
  shape = c("even", "bimodal", "concentrated_precise", "concentrated_noisy"),
  sigma_tau = 1
)
```

## Arguments

- J:

  Number of sites, an integer of at least 1.

- I:

  Informativeness of the profile, a number strictly between 0 and 1.

- V:

  Logit variance of the profile, a number of at least 0. With `J = 1` it
  must be 0.

- shape:

  Shape of the profile: one of `"even"` (the default), `"bimodal"`,
  `"concentrated_precise"` or `"concentrated_noisy"`. See Details.

- sigma_tau:

  Between-site standard deviation at which the profile has the
  informativeness `I`, a positive number.

## Value

A numeric vector of `J` sampling variances.

## Details

The variances are \\v_j = \tilde v\\ e^{x_j}\\ with \\\tilde v =
\sigma\_\tau^2 (1 - I) / I\\ and deviations \\x_j\\ of mean 0 and mean
square \\V\\, so that the logits of the reliabilities have mean
\\\mathrm{logit}(I)\\ and variance \\V\\. The shapes of the deviations:

- `"even"`: equally spaced, a ladder of variances that is even on the
  log scale.

- `"bimodal"`: two groups of \\\lceil J/2 \rceil\\ precise and \\\lfloor
  J/2 \rfloor\\ noisy sites, which are not symmetric when \\J\\ is odd.

- `"concentrated_precise"`: \\k = \lceil J/5 \rceil\\ precise sites at
  \\-\sqrt{V(J-k)/k}\\ and the other \\J - k\\ at \\\sqrt{Vk/(J-k)}\\,
  so a fifth of the sites have reliabilities far above the rest.

- `"concentrated_noisy"`: the mirror image, with \\k\\ noisy sites.

At `V = 0` all four shapes give equal variances. The sites come in the
order of their deviations as listed: from precise to noisy, except for
`"concentrated_noisy"`, whose noisy sites come first.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
v <- shape_profile(J = 40, I = 0.3, V = 1, shape = "concentrated_precise")
S <- site_reliability(se2 = v, sigma_tau = 1)
kl_centroid(S)
#> [1] 0.3
logit_moments(S)
#>    V rho3 
#>  1.0  1.7 
```
