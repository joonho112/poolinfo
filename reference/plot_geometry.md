# Geometry of the Kullback-Leibler centroid

Shows the geometry behind the informativeness as a centroid of the
reliabilities: (a) the Bernoulli negative entropy, a tangent line at
each of a few sites, and the vertical gap between the curve and each
tangent at \\p = I\\, which is that site's term of the dispersion \\H\\;
(b) the same gaps as bars. Over all the sites, the mean of the gaps is
\\H\\.

## Usage

``` r
plot_geometry(x, sites = NULL, ...)
```

## Arguments

- x:

  An object of class `reliability_profile`, `reliability_posterior`,
  `reliability_draws` or `reliability_design`, or site reliabilities as
  a numeric vector.

- sites:

  Positions or labels of the sites to show, at most 6. By default four
  sites spread over the range of the reliabilities.

- ...:

  Must be empty.

## Value

A patchwork object with two panels.

## Details

Site \\j\\'s term of \\H\\,
\\\mathrm{KL}\\\mathrm{Ber}(I)\\\mathrm{Ber}(S_j)\\\\, is the Bregman
divergence of \\\phi(p) = p \log p + (1 - p)\log(1 - p)\\: the height of
\\\phi\\ at \\I\\ above its tangent at \\S_j\\ (Lee, 2026, Section 4).
For a posterior fit the profile is taken at the posterior median of
`sigma_tau`. Panel (b) shows these sites only; their mean is not \\H\\
unless all the sites are shown.

## Examples

``` r
plot_geometry(c(0.2, 0.4, 0.6, 0.85))
```
