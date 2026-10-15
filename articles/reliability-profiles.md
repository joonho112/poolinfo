# Reliability profiles and their summaries

Two trials can have the same informativeness and still differ in what
partial pooling does to their sites. In one, every site is about equally
reliable. In another, a few precise sites stand next to many noisy ones.
In this vignette we take the between-site standard deviation
$`\sigma_\tau`$ as known and study the reliabilities themselves. We look
at how the informativeness averages them, what the dispersion and the
concentration count add, and how all four summaries change with
$`\sigma_\tau`$, using the geometry of Lee ([2026](#ref-lee2026))
without deriving it.

``` r

library(poolinfo)
```

Site $`j`$’s reliability is
$`S_j = \sigma_\tau^2 / (\sigma_\tau^2 + v_j)`$, where $`v_j`$ is the
sampling variance of the site’s estimate, and a trial’s reliability
profile is the set of its $`J`$ reliabilities. See
[`vignette("poolinfo")`](https://joonho112.github.io/poolinfo/articles/poolinfo.md)
for why $`S_j`$ is both a reliability and the shrinkage weight of the
normal model.

## Five profiles with the same informativeness

Our examples come from the data set `example_profiles`, with the
sampling variances of five artificial trials of 40 sites each and the
reliabilities that they give at $`\sigma_\tau = 0.2`$. All five have the
informativeness $`I = 0.3`$. They differ in their logit variance $`V`$,
the variance of the logits $`\mathrm{logit}\, S_j`$ with divisor $`J`$.
Because $`\mathrm{logit}\, S_j = \log\sigma_\tau^2 - \log v_j`$, $`V`$
is also the variance of the log sampling variances, and it does not
change with $`\sigma_\tau`$. The five profiles are

- `homogeneous`, in which every site has the same sampling variance, so
  that $`V = 0`$;
- `even`, with log variances spread evenly;
- `bimodal`, with two groups of 20 sites, one precise and one noisy;
- `concentrated_precise`, with 8 precise sites and 32 noisy ones;
- `concentrated_noisy`, with 32 precise sites and 8 noisy ones.

Apart from `homogeneous`, all have $`V = 1.5`$. R’s
[`var()`](https://rdrr.io/r/stats/cor.html) divides by $`J - 1`$ instead
and gives 1.54 for them.

We make a profile of each trial with
[`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md),
from its sampling variances and $`\sigma_\tau`$, and collect the five
profiles with
[`compare_profiles()`](https://joonho112.github.io/poolinfo/reference/compare_profiles.md).

``` r

by_shape <- split(example_profiles, example_profiles$profile)
profiles <- lapply(by_shape, function(d) {
  reliability_profile(se2 = d$se2, sigma_tau = 0.2)
})
cmp <- compare_profiles(.list = profiles)
cmp
#> Comparison of 5 analyses
#> Each analysis is conditional on its own standard errors.
#> 
#>   analysis              J   informativeness I  mean reliability
#>   homogeneous           40  0.300              0.300
#>   even                  40  0.300              0.347
#>   bimodal               40  0.300              0.353
#>   concentrated_precise  40  0.300              0.317
#>   concentrated_noisy    40  0.300              0.360
#> 
#>   analysis              dispersion H  concentration count J_eff
#>   homogeneous           0.0000        40.0
#>   even                  0.1486        27.6
#>   bimodal               0.1524        27.3
#>   concentrated_precise  0.1676        24.1
#>   concentrated_noisy    0.1166        33.3
```

In every row $`I = 0.300`$, while the mean reliability runs from 0.300
in the homogeneous profile to 0.360. Although the four unequal profiles
have the same logit variance, their dispersions range from 0.117 to
0.168 and their concentration counts from 24.1 to 33.3, against 0 and 40
for the homogeneous profile. Of the four, the profile with 8 noisy sites
has the smallest dispersion and yet the largest mean reliability, 0.360.
Its mirror image, with 8 precise sites, combines the largest dispersion
with the smallest mean reliability among them, 0.317. So a larger
dispersion need not come with a larger gap between the mean reliability
and $`I`$([Lee 2026](#ref-lee2026)). In every column the even and
bimodal profiles are close to each other.

Sorting the reliabilities of the two concentrated profiles shows where
these differences come from. We put the two in one comparison, to use
again below, and plot their sorted reliabilities in one panel.

``` r

two <- compare_profiles(
  `8 precise` = profiles$concentrated_precise,
  `8 noisy` = profiles$concentrated_noisy
)
plot(two, type = "profile")
```

![The sorted site reliabilities of the profiles with 8 precise and with
8 noisy sites, in one panel. In the first, 32 sites lie below the
informativeness and 8 far above it; in the second, 8 lie far below it
and 32 above it. A horizontal line marks the informativeness of both,
0.3.](reliability-profiles_files/figure-html/two-shapes-1.png)

With 8 precise sites, the other 32 sit together at 0.19, below $`I`$,
and the precise ones at 0.83. In the profile with 8 noisy sites it is
the other way around, with 32 sites at 0.44, above $`I`$, and 8 at 0.04.
In each profile the few sites lie far from $`I`$ on the logit scale and
the many lie close to it on the other side, so that the logits average
to $`\mathrm{logit}\, I`$.

## Three averages of the reliabilities

On the logit scale, the informativeness is the plain average of the
reliabilities,

``` math
\mathrm{logit}\, I = \frac{1}{J}\sum_{j=1}^J \mathrm{logit}\, S_j,
\qquad
I = \frac{\sigma_\tau^2}{\sigma_\tau^2 + \mathrm{GM}(v)},
```

where $`\mathrm{GM}(v)`$ is the geometric mean of the sampling
variances. Because the logits are $`\log\sigma_\tau^2 - \log v_j`$,
averaging them averages the log variances, which gives the second form.

Lee ([2026](#ref-lee2026), Theorem 1) shows that $`I`$ is also the value
$`p`$ that minimizes the average Kullback-Leibler divergence

``` math
\frac{1}{J}\sum_{j=1}^J \mathrm{KL}\{\mathrm{Ber}(p) \,\|\, \mathrm{Ber}(S_j)\}.
```

We call the minimizer of such an average a centroid of the
reliabilities. This divergence is one of the Bregman divergences, and
two others give the other familiar means. With $`p`$ in the first
argument, as above, squared error gives the arithmetic mean $`\bar S`$,
and the Bregman divergence of the Burg potential $`-\log p`$ gives the
harmonic mean, our harmonic centroid. With $`p`$ in the second argument,
all three divergences give $`\bar S`$.
[`reliability_centroids()`](https://joonho112.github.io/poolinfo/reference/reliability_centroids.md)
returns the three first-argument centroids, and we apply it to each
profile.

``` r

centroids <- sapply(by_shape, function(d) reliability_centroids(d$reliability))
round(t(centroids), 3)
#>                      harmonic   I  mean
#> homogeneous             0.300 0.3 0.300
#> even                    0.181 0.3 0.347
#> bimodal                 0.188 0.3 0.353
#> concentrated_precise    0.223 0.3 0.317
#> concentrated_noisy      0.135 0.3 0.360
```

In the four unequal profiles the harmonic centroid, from 0.135 to 0.223,
lies below $`I`$. It is never larger than the other two, in any profile,
and the three are equal only when all the reliabilities are, as in the
homogeneous profile. Between $`I`$ and $`\bar S`$, however, there is no
general order. Here $`\bar S > I`$ in all four unequal profiles, the
sign of $`1 - 2I`$. A profile close enough to constant must have that
sign when $`I \ne 1/2`$, but a profile far from constant can have the
other. Lee ([2026](#ref-lee2026)) gives one with 25 sites and
$`I = 0.3`$, where one site is far more precise than the other 24, which
are all alike. We build it from its centered logits
$`\Delta_j = \mathrm{logit}\, S_j - \mathrm{logit}\, I`$, with mean 0
and variance $`V = 2`$.

``` r

delta <- c(sqrt(48), rep(-sqrt(1 / 12), 24))
S <- plogis(qlogis(0.3) + delta)
reliability_centroids(S)
#>  harmonic         I      mean 
#> 0.2506438 0.3000000 0.2732468
```

Its mean reliability, 0.273, is below $`I`$ although $`I < 1/2`$. In the
next section we give a condition that guarantees the sign of $`1 - 2I`$.

## The dispersion and the gap between the averages

The dispersion $`H`$ is the smallest value of the average divergence
above, reached at $`p = I`$,

``` math
H = \frac{1}{J}\sum_{j=1}^J \mathrm{KL}\{\mathrm{Ber}(I) \,\|\, \mathrm{Ber}(S_j)\},
```

the average divergence of $`I`$ from the site reliabilities. It is 0
only for a constant profile.

The sign of the gap $`\bar S - I`$ matters for the site effects. If
$`\sigma_\tau`$ and the mean effect were known, the normal model would
give the site effects an average posterior variance of
$`\sigma_\tau^2(1 - \bar S)`$([Lee 2026](#ref-lee2026), Proposition 2).
Equally reliable sites with the same $`I`$ would give
$`\sigma_\tau^2(1 - I)`$, so a positive gap means a smaller average
posterior variance than theirs, and a negative gap a larger one.

Both $`H`$ and the gap have leading terms proportional to the logit
variance $`V`$. With $`\rho_3 = J^{-1}\sum_j |\Delta_j|^3`$, the third
absolute moment of the centered logits, Proposition 2 of Lee
([2026](#ref-lee2026)) gives

``` math
\begin{aligned}
H &= \tfrac12 I(1-I)V + R_H, &
|R_H| &\le \tfrac{\sqrt3}{108}\,\rho_3,\\
\bar S - I &= \tfrac12 I(1-I)(1-2I)V + R_S, &
|R_S| &\le \tfrac{1}{48}\,\rho_3.
\end{aligned}
```

When the gap’s leading term is larger than its bound in absolute value,
the gap has the sign of $`1 - 2I`$.
[`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md)
gives, for $`H`$ and for the gap, the exact value, the leading term and
the bound, and a logical column, `holds`, that is `TRUE` when the
leading term is the larger. We look first at the bimodal profile.

``` r

local_terms(by_shape$bimodal$reliability)
#>   quantity      exact leading      bound holds
#> 1        H 0.15241610  0.1575 0.02946278  TRUE
#> 2      gap 0.05255058  0.0630 0.03827328  TRUE
```

For this profile the gap’s leading term, 0.063, is larger than its
bound, 0.038, so the gap is positive, as $`1 - 2I`$ is, and its exact
value is 0.053. Since $`H`$ is never negative, its bound measures
instead how close the leading term is. Here $`H = 0.152`$, within 0.030
of 0.158.

With 8 precise sites the leading terms are the same, because $`I`$ and
$`V`$ are, but a few of the centered logits are large, and $`\rho_3`$
rises from 1.84 to 3.12.

``` r

local_terms(by_shape$concentrated_precise$reliability)
#>   quantity      exact leading      bound holds
#> 1        H 0.16758423  0.1575 0.05008673  TRUE
#> 2      gap 0.01727911  0.0630 0.06506457 FALSE
```

Now the gap’s bound, 0.065, is larger than its leading term, 0.063, and
`holds` is `FALSE`, so the leading term alone does not fix the sign. The
gap is still positive, but at 0.017 it is only 27% of the leading term.
Its remainder, -0.046, cancels most of the leading term and comes from
the 8 precise sites, whose logits lie far above the rest. In the profile
of 25 sites the bound is 0.278, more than three times the leading term,
and there the remainder reverses the sign.

Site $`j`$’s term of $`H`$,
$`\mathrm{KL}\{\mathrm{Ber}(I) \,\|\, \mathrm{Ber}(S_j)\}`$, is the
divergence of $`I`$ from that site’s reliability.
[`kl_contributions()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md)
returns the terms, and
[`plot_contributions()`](https://joonho112.github.io/poolinfo/reference/plot_contributions.md)
plots them in increasing order, here for the profile with 8 precise
sites.

``` r

plot_contributions(profiles$concentrated_precise)
```

![Each site's term of the dispersion H for the profile with 8 precise
sites, in increasing order, with H as a dashed line. The 32 other sites
have small terms and the 8 precise sites much larger
ones.](reliability-profiles_files/figure-html/contributions-1.png)

Each of the 8 precise sites has a term of 0.694, against 0.036 for each
of the other 32, so these 8 sites account for 83% of $`H`$.

## Changing the between-site standard deviation

The sampling variances fix the centered logits,
$`\Delta_j = \log \mathrm{GM}(v) - \log v_j`$, and $`\sigma_\tau`$ moves
all the logits together. As $`\sigma_\tau`$ grows from 0, $`I`$ grows
from 0 to 1, and along this path the derivative of $`H`$ with respect to
the logit of $`I`$ is the gap $`\bar S - I`$([Lee 2026](#ref-lee2026),
Proposition 2). So $`H`$ rises from 0 while the mean reliability is
above $`I`$, reaches its largest value where the gap is 0, and falls
back to 0. To see how far the summaries move, we compute the profile
with 8 precise sites at four values of $`\sigma_\tau`$ with
[`update()`](https://rdrr.io/r/stats/update.html).

``` r

precise <- profiles$concentrated_precise
sigmas <- c(0.15, 0.2, 0.25, 0.3)
along <- sapply(sigmas, function(s) {
  update(precise, sigma_tau = s)$stats[c("I", "mean", "gap", "H", "J_eff")]
})
round(cbind(sigma_tau = sigmas, t(along)), 3)
#>      sigma_tau     I  mean    gap     H  J_eff
#> [1,]      0.15 0.194 0.240  0.045 0.149 19.298
#> [2,]      0.20 0.300 0.317  0.017 0.168 24.114
#> [3,]      0.25 0.401 0.390 -0.011 0.169 28.505
#> [4,]      0.30 0.491 0.458 -0.033 0.161 31.959
```

Over this range $`H`$ stays between 0.149 and 0.169, while the gap
changes sign, from 0.045 to -0.033, and the concentration count rises
from 19.3 to 32.0. At $`\sigma_\tau = 0.25`$ the mean reliability,
0.390, is below $`I = 0.401`$ although $`I < 1/2`$, as in the profile of
25 sites.

[`profile_path()`](https://joonho112.github.io/poolinfo/reference/profile_path.md)
follows the whole path with the centered logits held fixed, and for a
profile it also gives $`\sigma_\tau`$ at each point. On a fine grid we
find where $`H`$ is largest.

``` r

path <- profile_path(precise, n = 2001)
path[which.max(path$H), c("sigma_tau", "I", "H", "gap")]
#>     sigma_tau         I         H          gap
#> 723 0.2296954 0.3611389 0.1700047 1.607891e-06
```

At the peak, $`H = 0.170`$, $`\sigma_\tau \approx 0.23`$ and
$`I \approx 0.36`$, and the gap is 0 up to the spacing of the grid. So
the peak is not at $`I = 1/2`$, where the gap’s leading term changes
sign.
[`plot_gap_path()`](https://joonho112.github.io/poolinfo/reference/plot_gap_path.md)
plots $`H`$ and the gap along the path, with their leading terms.

``` r

plot_gap_path(precise)
```

![Two panels along the informativeness from 0 to 1 for the profile with
8 precise sites. The dispersion H peaks before I = 1/2, at the point
where the gap between the mean reliability and I crosses 0. The dashed
leading terms peak and cross 0 at I =
1/2.](reliability-profiles_files/figure-html/path-1.png)

In the left panel $`H`$ rises steeply from 0, peaks at the cross and
returns to 0 at $`I = 1`$, while its leading term (dashed) peaks at
$`I = 1/2`$. In the right panel the gap crosses 0 at the same $`I`$,
about 0.36, and stays negative up to $`I = 1`$. For every $`I`$ strictly
between this peak and 1/2, this profile has $`\bar S < I`$, against the
sign of its leading term. The point on each curve is the profile at
$`\sigma_\tau = 0.2`$, before the peak. A much larger $`\sigma_\tau`$
would bring all the reliabilities close to 1, where $`H`$ and the gap
return to 0. Along every path with unequal sampling variances the gap
changes sign once, because $`H`$ is strictly concave in $`I`$ when the
centered logits are fixed. With equal sampling variances both $`H`$ and
the gap are 0 throughout.
[`vignette("implementation")`](https://joonho112.github.io/poolinfo/articles/implementation.md)
gives the argument.

## The concentration count

With the proportions $`q_j = S_j / \sum_k S_k`$ of the total
reliability, the concentration count is the inverse of the sum of their
squares,

``` math
J_{\mathrm{eff}} = \frac{(\sum_j S_j)^2}{\sum_j S_j^2}
= \frac{1}{\sum_j q_j^2}.
```

Its largest value, $`J`$, belongs to profiles of equally reliable sites,
and it approaches 1, without reaching it, as one site’s proportion
approaches 1. It depends only on the proportions. A trial whose 40 sites
all have reliability 0.05 and one whose 40 sites all have 0.95 both have
$`J_{\mathrm{eff}} = 40`$.

We plot the cumulative proportions of the two concentrated profiles, the
most reliable sites first, with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
`type = "concentration"`, a call to
[`plot_concentration()`](https://joonho112.github.io/poolinfo/reference/plot_concentration.md).
Because the 8 precise sites are tied, we leave out the labels that would
name the most reliable sites.

``` r

plot(two, type = "concentration", label_top = 0)
```

![The cumulative proportion of the total reliability against the
fraction of sites, most reliable first, for the profiles with 8 precise
and with 8 noisy sites, with the diagonal of equal proportions. The
first curve climbs steeply over the first fifth of the sites; the second
climbs a little faster than the diagonal and flattens after four fifths
of the
sites.](reliability-profiles_files/figure-html/concentration-1.png)

With 8 precise sites the curve climbs steeply over the first fifth of
the sites, which carry 52% of the total reliability, and then rises
slowly. With 8 noisy sites it climbs evenly, a little faster than the
diagonal, until the 32 precise sites have contributed 98% of the total,
and the 8 noisy ones add the rest. The concentration counts printed in
the panel, 24.1 and 33.3, summarize the two curves. To see which sites
carry the reliability in your own analysis, you can give a posterior fit
to
[`plot_concentration()`](https://joonho112.github.io/poolinfo/reference/plot_concentration.md)
and
[`plot_contributions()`](https://joonho112.github.io/poolinfo/reference/plot_contributions.md)
as well.

## The geometry in two plots

Each site’s term of $`H`$ is a Bregman divergence. With
$`\phi(p) = p \log p + (1 - p)\log(1 - p)`$, the term
$`\mathrm{KL}\{\mathrm{Ber}(I) \,\|\, \mathrm{Ber}(S_j)\}`$ is the
height of $`\phi`$ at $`I`$ above the tangent to $`\phi`$ at $`S_j`$.
[`plot_geometry()`](https://joonho112.github.io/poolinfo/reference/plot_geometry.md)
shows the curve, the tangents at a few sites and these heights. We use
the even profile, whose 40 reliabilities all differ, and four of its
sites spread over the range.

``` r

even <- profiles$even
plot_geometry(even, sites = c(1, 14, 27, 40))
```

![Panel (a): the curve p log p + (1 - p) log(1 - p) with the tangents at
four sites of the even profile and, at p = I, the vertical gap between
the curve and each tangent. Panel (b): the four gaps as
bars.](reliability-profiles_files/figure-html/geometry-1.png)

In panel (a) the longest gap belongs to the most reliable of the four
sites, at 0.77, and the next longest to the least reliable, at 0.05.
These two sites are equally far from $`I`$ on the logit scale, but the
divergence is not symmetric, and the site above $`I`$ has the larger
term. Between them, the sites at 0.46 and 0.18 are closer to $`I`$ and
have short gaps. Panel (b) gives the four terms as bars, 0.503, 0.054,
0.045 and 0.317. Their mean is not $`H`$, the average over all 40 sites.

Theorem 2 of Lee ([2026](#ref-lee2026)) splits the average divergence of
a constant reliability $`p`$ from the profile into two parts,

``` math
\frac{1}{J}\sum_{j=1}^J \mathrm{KL}\{\mathrm{Ber}(p) \,\|\, \mathrm{Ber}(S_j)\}
= \mathrm{KL}\{\mathrm{Ber}(p) \,\|\, \mathrm{Ber}(I)\} + H,
```

so the average is smallest at $`p = I`$, where it equals $`H`$.
[`plot_decomposition()`](https://joonho112.github.io/poolinfo/reference/plot_decomposition.md)
plots the left side against $`p`$ and shows the split at one value of
$`p`$, here 0.75.

``` r

plot_decomposition(even, at = 0.75)
```

![The average divergence of a constant reliability p from the even
profile against p, with its minimum H at p = I, the mean reliability
marked just to its right, and the split at p = 0.75 into the divergence
of p from I plus
H.](reliability-profiles_files/figure-html/decomposition-1.png)

The curve’s minimum, at $`I = 0.3`$, is $`H = 0.149`$, the height of the
dashed line. At $`p = 0.75`$ the average divergence is 0.578, of which
0.430 is the divergence of 0.75 from $`I`$ and the rest is $`H`$. The
mean reliability, 0.347, sits just to the right of the minimum and a
little above it on the curve. It minimizes the average with $`p`$ in the
second argument, a different function of $`p`$.

## Making profiles of your own

You can build profiles like these with
[`shape_profile()`](https://joonho112.github.io/poolinfo/reference/shape_profile.md),
for a simulation or to see what the diagnostics of a planned trial might
look like. It returns sampling variances with a given number of sites,
informativeness, logit variance and shape at a given $`\sigma_\tau`$.
Here we ask for 20 sites in two groups, with $`I = 1/2`$ and $`V = 1`$.

``` r

v <- shape_profile(J = 20, I = 0.5, V = 1, shape = "bimodal",
                   sigma_tau = 0.25)
half <- reliability_profile(se2 = v, sigma_tau = 0.25)
half
#> Reliability profile of 20 sites at sigma_tau = 0.25 (supplied)
#> 
#>   informativeness I            0.500
#>   mean reliability             0.500
#>   dispersion H                 0.1201
#>   logit variance V             1.000
#>   concentration count J_eff    16.5 of 20
#> 
#> Site reliabilities range from 0.27 to 0.73.
```

Its reliabilities take two values, 0.27 and 0.73, whose logits lie
symmetrically around $`\mathrm{logit}\, I = 0`$. At $`I = 1/2`$ the
gap’s leading term is 0, and for this symmetric profile the mean
reliability equals $`I`$. Yet $`H = 0.120`$, close to its leading term
$`\tfrac12 I(1 - I)V = 0.125`$, so two equal centroids do not make the
reliabilities equal.

[`vignette("planning")`](https://joonho112.github.io/poolinfo/articles/planning.md)
starts instead from the site sizes, and
[`vignette("posterior-summaries")`](https://joonho112.github.io/poolinfo/articles/posterior-summaries.md)
returns to an estimated $`\sigma_\tau`$.

## References

Lee, JoonHo. 2026. “The Geometry of Informativeness: A Bregman-Centroid
Diagnostic for Partial Pooling in Multisite Trials.” Unpublished
manuscript.
