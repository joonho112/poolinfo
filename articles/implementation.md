# Computation and reproducibility

In this vignette we describe how poolinfo computes its posterior
summaries and what evidence there is that they are accurate. We begin
with the marginal posterior of the between-site standard deviation and
the conditions under which its moments exist, then describe the default
quadrature rule and the checks that come with every fit. Next come the
second engine and how the two engines compare, the tests’ tolerances,
and our reproduction of the numbers in Lee ([2026](#ref-lee2026)). We
end with the derivations behind the parts of the package that go beyond
the paper, and a table that maps the paper’s notation to the package’s.

``` r

library(poolinfo)
```

## Integrating out the mean

Collect the site estimates $`\hat\tau_1, \dots, \hat\tau_J`$ in a vector
$`y`$, and let $`X`$ be the $`J \times p`$ design matrix of the mean
model, with rows $`x_j^\top`$. For a common mean it is a column of ones,
and $`p = 1`$. Site effects are normal around $`x_j^\top\beta`$ with
standard deviation $`\sigma_\tau`$, and the estimates are normal around
the effects with known variances $`v_j`$. Write
$`w_j = 1/(\sigma_\tau^2 + v_j)`$ and
$`W = \mathrm{diag}(w_1, \dots, w_J)`$, let
$`b = (X^\top W X)^{-1} X^\top W y`$ be the weighted least-squares fit,
and let $`Q = (y - Xb)^\top W (y - Xb)`$ be its weighted residual sum of
squares. Under a flat prior on $`\beta`$, integrating $`\beta`$ out of
the likelihood leaves the posterior density of $`\sigma_\tau`$, up to a
constant,

``` math
g(\sigma_\tau) = \pi(\sigma_\tau)\, |X^\top W X|^{-1/2}
\prod_j (\sigma_\tau^2 + v_j)^{-1/2} \exp(-Q/2),
```

where $`\pi`$ is the prior density of $`\sigma_\tau`$. Apart from
$`\pi`$, $`g`$ is the restricted likelihood of a random-effects
meta-regression at $`\tau^2 = \sigma_\tau^2`$ in metafor’s notation
([Harville 1974](#ref-harville1974); [Viechtbauer
2010](#ref-viechtbauer2010)). Under the flat prior on $`\sigma_\tau`$,
then, the posterior density of $`\sigma_\tau`$ peaks at the REML
estimate of $`\sigma_\tau`$. For a common mean,
$`X^\top W X = \sum_j w_j`$, $`b`$ is the precision-weighted mean, and
$`g`$ is equation (D.2) of Lee ([2026](#ref-lee2026)). Every diagnostic
is a function of $`\sigma_\tau`$ and the $`v_j`$ alone, so each of its
posterior summaries is a one-dimensional integral against $`g`$.

## Tails and moments

For every $`\sigma_\tau \ge B > 0`$ we have three bounds, where
$`v_{\max}`$ is the largest $`v_j`$:

- $`\exp(-Q/2) \le 1`$,
- $`\prod_j (\sigma_\tau^2 + v_j)^{-1/2} \le \sigma_\tau^{-J}`$, and
- $`X^\top W X \ge X^\top X / (\sigma_\tau^2 + v_{\max})`$ in the matrix
  order.

Together they give

``` math
g(\sigma_\tau) \le |X^\top X|^{-1/2} (1 + v_{\max}/B^2)^{p/2}\,
\sigma_\tau^{p - J}\, \pi(\sigma_\tau), \qquad \sigma_\tau \ge B,
```

which for $`p = 1`$ is the tail bound of Appendix D.3 of the paper. In
the limit the bound is sharp. As $`\sigma_\tau`$ grows, $`Q`$ tends to 0
and the ratio of $`g(\sigma_\tau)`$ to
$`\sigma_\tau^{p - J} \pi(\sigma_\tau)`$ tends to $`|X^\top X|^{-1/2}`$,
so the tail of $`g`$ is that of
$`\sigma_\tau^{p - J} \pi(\sigma_\tau)`$.

Suppose the prior density falls off like $`\sigma_\tau^{-k}`$, with
$`k = 0`$ for the flat prior, $`k = \nu + 1`$ for a half-t prior with
$`\nu`$ degrees of freedom and $`k = 2`$ for the half-Cauchy. Then the
posterior mean of $`\sigma_\tau^q`$ is finite exactly when
$`J > q + p + 1 - k`$, and the case $`q = 0`$ is the condition for a
proper posterior. Under the half-normal and uniform priors every moment
exists. A flat prior thus needs $`J \ge p + 2`$ sites. A proper prior
gives a proper posterior with as few as $`p`$ sites. With $`J = p`$,
however, the fitted values reproduce the estimates whatever
$`\sigma_\tau`$ is, $`g`$ is a constant times $`\pi`$, and the posterior
is the prior. For that reason the package asks for at least $`p + 1`$
sites under a proper prior. It integrates only the moments that exist
and reports the others as not defined. A posterior mean of a mean-model
coefficient requires the first absolute moment of its marginal
distribution. Since its conditional standard deviation grows
proportionally to $`\sigma_\tau`$, this requires $`E(\sigma_\tau)`$ to
exist: at least $`p + 3`$ sites under the flat prior. Its standard
deviation requires $`E(\sigma_\tau^2)`$, hence $`p + 4`$ sites. Central
intervals remain available when these moments do not exist.

The upper end $`B`$ of the integrals is doubled until the bound on the
omitted tail of every moment that exists is within its limit, at most
400 times. A moment that exists can still need a larger $`B`$ than that
when its tail falls off very slowly: the mean of $`\sigma_\tau`$ with 2
sites and a half-t prior with 0.05 degrees of freedom has a tail bound
that falls like $`B^{-0.05}`$. Both engines then stop with an error that
names the moment. A prior with a lighter tail avoids it.

## Quadrature on panels

The default engine, `posterior_control(engine = "nodes")`, builds one
quadrature rule in $`\sigma_\tau`$ for each fit and uses it for every
summary of the fit. Above a lower end $`\sigma_{\mathrm{lo}}`$ the rule
is composite Gauss-Legendre quadrature on $`t = \log\sigma_\tau`$, with
the Jacobian $`d\sigma_\tau = \sigma_\tau\,dt`$ in its weights. At the
start no panel is wider than 0.5 on the log scale. We halve a panel
until its 16-point and 32-point rules agree, for the posterior mass and
for the second moment when it exists. They agree when they differ by
less than $`10^{-14}`$ of the total or $`64J\epsilon`$ of the panel’s
own integral, whichever is larger, where $`\epsilon`$ is the machine
epsilon. Below the second limit, the rounding error of $`\log g`$, a sum
of $`J`$ terms, no amount of halving brings the two rules closer. Each
panel of the final rule has 32 nodes.

Below $`\sigma_{\mathrm{lo}}`$ the rule is a single 32-point
Gauss-Legendre panel in $`\sigma_\tau`$ itself. Panels on the log scale
cannot reach 0, yet under every prior of
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md)
the posterior density of $`\sigma_\tau`$ is positive there, and roughly
$`g(0)\,\sigma_{\mathrm{lo}}`$ of the unnormalized mass lies between 0
and $`\sigma_{\mathrm{lo}}`$. We set $`\sigma_{\mathrm{lo}}`$ to
$`e^{-8}`$ times the smallest standard error, or times the scale of a
proper prior or the upper limit of a uniform prior if that is smaller.
Below it $`\sigma_\tau^2`$ is less than $`e^{-16}`$ times every $`v_j`$,
and neither the likelihood nor the prior changes much.

Under a uniform prior the upper end $`B`$ is the prior’s upper limit.
Otherwise $`B`$ starts well beyond the spread of the estimates and the
standard errors, and is doubled until the bounds of the previous
section, divided by the mass below $`B`$, meet two limits. The bound on
the mass must fall below $`10^{-12}`$, and the bound on
$`E(\sigma_\tau^q)`$ below $`10^{-10}\,\mathrm{GM}(v)^{q/2}`$, for each
of $`q = 1, 2, 4`$ whose moment exists. Dividing by the mass below $`B`$
rather than the whole makes the bounds conservative. With the factor
$`\mathrm{GM}(v)^{q/2}`$ the limits change with the units of the
estimates, so the rule does not depend on them. We compare the bounds
with their limits on the log scale, where no power of $`\sigma_\tau`$
overflows. Within a panel, the distribution function is the integral of
the Legendre series that matches the density at the panel’s 32 nodes,
and the quantiles are its roots.

We use the passing outcome of the credit recovery study as the example.
Each fit stores its rule as `fit$rule`.

``` r

fit <- reliability_posterior(estimate, se = se, data = ocrs)
rule <- fit$rule
data.frame(nodes = length(rule$sigma), panels = length(rule$edges) - 1,
           sigma_lo = rule$sigma_lo, B = rule$upper)
#>   nodes panels     sigma_lo        B
#> 1   800     24 2.596037e-05 3.630941
```

The rule has 800 nodes, 32 in each of the 24 panels above
$`\sigma_{\mathrm{lo}}`$, here 2.6e-05, and 32 below it. Its upper end,
$`B = 3.63`$, is 28 times $`\sqrt{\mathrm{GM}(v)} = 0.128`$ and 14 times
the 97.5% quantile of $`\sigma_\tau`$.

## Numerical checks

Every fit carries checks of its integration, which
[`check_numerics()`](https://joonho112.github.io/poolinfo/reference/check_numerics.md)
lists.

``` r

check_numerics(fit)[, c("check", "value", "limit", "status")]
#>                                                             check        value
#> 1     normalizing constant, rule against halved panels (relative) 0.000000e+00
#> 2            E(I) and E(H), rule against halved panels (relative) 1.654837e-16
#> 3           normalizing constant, lower end moved down (relative) 0.000000e+00
#> 4                  E(I) and E(H), lower end moved down (relative) 0.000000e+00
#> 5                               mass beyond the upper end (bound) 9.678885e-32
#> 6 E(sigma_tau^1) beyond the upper end, over sqrt(GM(v))^1 (bound) 2.844392e-30
#> 7 E(sigma_tau^2) beyond the upper end, over sqrt(GM(v))^2 (bound) 8.368938e-29
#> 8 E(sigma_tau^4) beyond the upper end, over sqrt(GM(v))^4 (bound) 7.273389e-26
#> 9            turning points of H at two resolutions (logit scale) 6.938894e-17
#>   limit status
#> 1 1e-10     ok
#> 2 1e-10     ok
#> 3 1e-10     ok
#> 4 1e-10     ok
#> 5 1e-12     ok
#> 6 1e-10     ok
#> 7 1e-10     ok
#> 8 1e-10     ok
#> 9 1e-07     ok
```

The first four rows compare the rule with two others. Halving every
panel tests whether the panels are fine enough, and moving the lower end
down to $`\sigma_{\mathrm{lo}}/2`$ tests the treatment of the mass near
0. Each compares the normalizing constants and the posterior means of
$`I`$ and $`H`$, as relative differences, with the limit $`10^{-10}`$,
or $`64J\epsilon`$ if that is larger. The next rows bound the mass
beyond $`B`$ and, for the orders that exist, the part of
$`E(\sigma_\tau^q)`$ beyond $`B`$, divided by
$`\sqrt{\mathrm{GM}(v)}^q`$. The last row compares $`H`$’s turning
point, on the logit scale, as found by two scans of different
resolution. For this fit the largest relative difference between rules
is 1.7e-16, and the tail bounds are at most 7.3e-26. Since $`B`$ is
doubled only until the tail bounds pass, in other fits those rows can
come close to their limits. A value outside its limit makes
[`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
give a warning that points to
[`check_numerics()`](https://joonho112.github.io/poolinfo/reference/check_numerics.md).

## The adaptive engine

`posterior_control(engine = "adaptive")` computes each posterior mean,
probability and quantile as a separate integral with
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html), as in
the computation of Lee ([2026](#ref-lee2026)). On a finite range
[`integrate()`](https://rdrr.io/r/stats/integrate.html) uses the
adaptive Gauss-Kronrod routine of QUADPACK ([Piessens et al.
1983](#ref-piessens1983)), which starts from a single 21-point rule over
the whole range. When the integrand is concentrated in a small part of a
long range, that first rule can miss it, and
[`integrate()`](https://rdrr.io/r/stats/integrate.html) returns 0 with
an error estimate of 0:

``` r

integrate(dnorm, 0, 1e4, mean = 60, sd = 0.5)
#> 0 with absolute error < 0
```

Its true value is 1. Splitting the range so that the peak falls in a
piece of moderate length avoids the failure. So the adaptive engine
splits the range at $`\sqrt{\mathrm{GM}(v)}`$, as the paper did, and
also at 9 points spread evenly on the log scale over the stretch where
the density of $`\log\sigma_\tau`$ exceeds $`e^{-70}`$ times its largest
value. Besides these extra splits, it differs from the paper’s
computation in two ways. Its integrals end at the node rule’s upper end
$`B`$, while the paper’s run to infinity. And it takes the posterior
standard deviations, and the intervals of the site effects and of the
coefficients, from a node rule. Each integrand is scaled before it is
integrated so that
[`integrate()`](https://rdrr.io/r/stats/integrate.html)’s absolute
tolerance acts as a relative one. For the eight schools, the tests check
that multiplying the estimates and the standard errors by factors from
$`10^{-100}`$ to $`10^{100}`$ changes the posterior means of the
diagnostics by less than $`10^{-9}`$, relative, under both engines.

The two engines use the same kernel $`g`$, the same upper end $`B`$ with
the tail bounds that fix it, and the same formulas for the diagnostics.
The adaptive engine also takes its standard deviations from a node rule.
So their agreement checks the integration over $`[0, B]`$, by panels and
Legendre series on one side and by adaptive subdivision on the other.
The kernel and the formulas are checked against metafor’s restricted
likelihood, against brute-force integrals, and against reference values
computed for Lee ([2026](#ref-lee2026)) with separate code that
integrates to infinity. The node rule is the default because it is
faster: a single rule serves every summary of a fit, while the adaptive
engine integrates each summary, and each site, separately.

## Tolerances

Below are the tolerances of the tests, kept in
`tests/testthat/helper-tolerances.R`. The first row is for identities of
the profile algebra, which are true in exact arithmetic. The others
govern two comparisons, with full-precision reference values computed
for Lee ([2026](#ref-lee2026)) and between the two engines. In the
fourth row, $`x`$ is the value being compared.

| Quantity | Tolerance |
|:---|:---|
| identities of the profile algebra | $`10^{-12}`$ |
| $`I`$, $`\bar S`$, $`H`$, $`S_j`$ and probabilities | $`10^{-9}`$, absolute |
| $`J_{\mathrm{eff}}`$ | $`10^{-9}`$, relative |
| $`\sigma_\tau`$, the mean effect, site means and standard deviations | $`10^{-9}\max\{\sqrt{\mathrm{GM}(v)}, \lvert x\rvert\}`$ |
| $`E(\sigma_\tau^2)`$ and its quantiles | $`10^{-8}`$, relative |
| $`I`$ at $`E(\sigma_\tau^2)`$ | $`10^{-8}`$, absolute |
| normalizing constant | $`10^{-10}`$, relative |
| turning point of $`H`$ (logit scale) | $`10^{-7}`$ |

For $`E(\sigma_\tau^2)`$ the tolerance is looser, because the posterior
tail of $`\sigma_\tau^2`$ is heavy. With 5 sites under the flat prior,
its variance does not exist. Brute-force integrals are less accurate
than either engine, and the tests that compare with them use looser
tolerances:

- $`10^{-8}`$ for posterior means under proper priors, and for the means
  and standard deviations of the site effects
- $`10^{-7}`$, relative, for $`E(\sigma_\tau)`$ when the prior’s scale
  is far below the standard errors
- $`10^{-6}`$ for the conditional moments of the site effects under a
  mean model, against a joint normal computation
- $`10^{-4}`$ for the kernel of a mean model, against a two-dimensional
  integral over $`\beta`$

## Reproducing the paper

The data set `ocrs` has the inputs of the credit recovery study’s 6
analyses in Lee ([2026](#ref-lee2026)). They are course passing in all
32 sites, in the 13 Algebra 1 sites and in the 19 English 9 sites,
passing with covariates, credits earned, and passing as a log odds
ratio. The first five are the rows of the paper’s Table 5 and Table G.1,
and the last is the log-odds row of Table G.2. We fit each with both
engines.

``` r

algebra <- ocrs$subject == "Algebra 1"
analyses <- list(
  Pooled = fit,
  `Algebra 1` = reliability_posterior(estimate, se = se,
                                      data = ocrs[algebra, ]),
  `English 9` = reliability_posterior(estimate, se = se,
                                      data = ocrs[!algebra, ]),
  Adjusted = reliability_posterior(estimate_adjusted, se = se_adjusted,
                                   data = ocrs),
  Credits = reliability_posterior(estimate_credits, se = se_credits,
                                  data = ocrs),
  `Log odds` = reliability_posterior(log_odds_ratio, se = se_log_odds_ratio,
                                     data = ocrs)
)
adaptive <- lapply(analyses, update,
                   control = posterior_control(engine = "adaptive"))
```

We tabulate the node rule’s posterior means and interval ends to four
decimals.

``` r

tab <- as.data.frame(compare_profiles(.list = analyses))
cols <- c("I", "I_lower", "I_upper", "mean", "H", "H_lower", "H_upper",
          "J_eff")
cbind(tab["analysis"], round(tab[cols], 4))
#>    analysis      I I_lower I_upper   mean      H H_lower H_upper   J_eff
#> 1    Pooled 0.6709  0.4974  0.8122 0.6635 0.0253  0.0185  0.0289 31.1704
#> 2 Algebra 1 0.6293  0.2939  0.8654 0.6253 0.0164  0.0086  0.0196 12.6976
#> 3 English 9 0.7074  0.4932  0.8678 0.6983 0.0284  0.0176  0.0342 18.4945
#> 4  Adjusted 0.6522  0.4736  0.8003 0.6444 0.0279  0.0205  0.0317 31.0343
#> 5   Credits 0.3181  0.0220  0.5972 0.3271 0.0285  0.0031  0.0376 28.7982
#> 6  Log odds 0.5302  0.2532  0.7434 0.5296 0.0346  0.0263  0.0367 30.1660
```

Rounded as the paper rounds them, the first five rows are the entries of
its Table 5, and the last row gives $`I`$ and $`H`$ in the log-odds row
of its Table G.2. The other numbers in the paper’s Tables G.1 and G.2
for these analyses, such as the posterior mean of $`\sigma_\tau`$, $`I`$
at $`E(\sigma_\tau^2)`$ and Pr(I \> 1/2), round to the paper’s entries
in the same way. Pooled passing, for example, has $`I = 0.6709`$ with
the interval from 0.4974 to 0.8122, which Table 5 prints as .67 \[.50,
.81\].

For each analysis we then find the largest difference between the two
engines over the posterior means and the ends of the intervals, in four
groups:

- the bounded diagnostics and the probability that $`I`$ exceeds 1/2
- $`J_{\mathrm{eff}}`$, relative to its value
- $`\sigma_\tau`$ and the mean effect, in units of
  $`\sqrt{\mathrm{GM}(v)}`$
- $`\sigma_\tau^2`$, relative to its value

``` r

differences <- function(a, b) {
  pa <- profile_summary(a)
  pb <- profile_summary(b)
  at <- function(p, q) unlist(p[p$quantity == q, c("mean", "lower", "upper")])
  gap <- function(q) abs(at(pa, q) - at(pb, q))
  rel <- function(q) abs(at(pa, q) / at(pb, q) - 1)
  c(bounded = max(gap("I"), gap("mean"), gap("H"),
                  abs(prob_above(a) - prob_above(b))),
    J_eff = max(rel("J_eff")),
    effect_scale = max(gap("sigma_tau"), gap("effect_mean")) /
      exp(mean(log(site_table(a)$se))),
    sigma_tau2 = max(rel("sigma_tau2")))
}
engines <- t(mapply(differences, analyses, adaptive))
signif(engines, 2)
#>           bounded   J_eff effect_scale sigma_tau2
#> Pooled    8.9e-16 2.2e-16      5.2e-15    5.1e-15
#> Algebra 1 1.1e-15 2.2e-16      8.3e-15    6.6e-15
#> English 9 3.3e-14 6.3e-15      2.5e-13    2.0e-13
#> Adjusted  2.3e-15 4.4e-16      1.1e-14    1.1e-14
#> Credits   6.2e-15 2.7e-15      1.9e-14    2.5e-13
#> Log odds  6.1e-16 2.2e-16      9.9e-16    3.1e-15
```

Across the 6 analyses the engines differ by at most 3.4e-14 on the
bounded diagnostics, 6.4e-15 relative on $`J_{\mathrm{eff}}`$, 2.6e-13
of $`\sqrt{\mathrm{GM}(v)}`$ on $`\sigma_\tau`$ and the mean effect, and
2.6e-13 relative on $`\sigma_\tau^2`$, at least 3 orders of magnitude
inside the tolerances.

The package’s tests make the same fits and compare them, under both
engines, with full-precision reference values computed for Lee
([2026](#ref-lee2026)) and stored in
`tests/testthat/cases/ocrs-reference.rds`. For each analysis they
compare 24 of the 28 stored summaries, including every number for these
analyses that the paper’s tables print, within the tolerances above. In
the same file, `tests/testthat/test-posterior-reference.R`, they do this
for 6 simulated data sets whose reference values are stored in
`tests/testthat/cases/boundary-cases.rds`. Slower tests in
`tests/testthat/test-engines-slow.R` compare the engines on the 6
analyses and on 500 random data sets of 3 to 200 sites. In these data
sets the sampling variances spread over several orders of magnitude, and
about a tenth of them have an outlying estimate. The slow tests also run
the node rule alone on 400 extreme data sets, with up to 10,000 sites
and ratios of sampling variances up to $`10^8`$, and require finite,
ordered summaries and every check within its limit. These tests take
several minutes and run only when the environment variable
`POOLINFO_SLOW` is `"true"`. From the package’s source directory, this
command runs them:

``` sh
POOLINFO_SLOW=true Rscript -e 'devtools::test(filter = "engines-slow")'
```

## Beyond the paper

The paper’s computations are for a common mean and the flat prior. The
package goes beyond them in five places, and for each we give the
derivation and the tests that check it.

### Mean models

With $`\beta`$ integrated out, site $`j`$’s effect given $`\sigma_\tau`$
has posterior mean $`S_j\hat\tau_j + (1 - S_j)\,x_j^\top b`$ and
variance
$`\sigma_\tau^2(1 - S_j) + (1 - S_j)^2\, x_j^\top (X^\top W X)^{-1} x_j`$,
where $`x_j^\top (X^\top W X)^{-1} x_j = h_{jj}/w_j`$ and $`h_{jj}`$ is
the $`j`$th diagonal element of $`W^{1/2}X`$’s hat matrix. The tests
compare $`g`$ with metafor’s restricted likelihood for $`p = 1, 2, 3`$
and with a two-dimensional brute-force integral over $`\beta`$. They
also compare the conditional moments of the site effects with a joint
normal computation in which $`\beta`$ has a very diffuse normal prior.

### Proper priors

The half-normal, half-t, half-Cauchy and uniform priors of
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md)
enter $`g`$ through $`\pi(\sigma_\tau)`$. The tail bound above contains
$`\pi(\sigma_\tau)`$, and the integral of the bound times
$`\sigma_\tau^q`$ beyond $`B`$ has a closed form, or a closed-form
bound, for each prior. It is a power of $`B`$ for the flat and half-t
priors, a normal tail probability or an upper incomplete gamma function
for the half-normal, and 0 for the uniform prior, whose upper limit is
$`B`$. We test the posteriors under each prior against brute-force
integrals, and each bound against the integral that it bounds.

### Intervals for H

With the $`v_j`$ fixed,
$`\mathrm{logit}\, S_j = \log\sigma_\tau^2 - \log v_j`$, so the centered
logits $`\Delta_j = \log\{\mathrm{GM}(v)/v_j\}`$ do not depend on
$`\sigma_\tau`$, which acts only through
$`m = \mathrm{logit}\, I = \log\{\sigma_\tau^2/\mathrm{GM}(v)\}`$. Along
this path $`dH/dm = \bar S - I`$([Lee 2026](#ref-lee2026), Proposition
2(i)), and differentiating once more gives

``` math
\frac{d^2 H}{dI^2} = -\frac{J^{-1}\sum_j (S_j - I)^2}{\{I(1 - I)\}^2},
```

which is negative unless all the $`v_j`$ are equal, and then $`H`$ and
the gap are 0 for every $`\sigma_\tau`$. With unequal sampling
variances, $`H`$ is therefore strictly concave in $`I`$ and has a single
peak, where $`\bar S = I`$. Lee ([2026](#ref-lee2026)) locates the peak
by scans at two resolutions; since there is only one, the package finds
the single root of $`\bar S - I`$ and keeps a scan as a check. For
$`0 < h < H_{\max}`$ there are two values $`\sigma_1(h) < \sigma_2(h)`$
of $`\sigma_\tau`$ at which $`H = h`$, and $`H \le h`$ if and only if
$`\sigma_\tau \le \sigma_1(h)`$ or $`\sigma_\tau \ge \sigma_2(h)`$. With
$`F`$ the posterior distribution function of $`\sigma_\tau`$,

``` math
P(H \le h) = F\{\sigma_1(h)\} + 1 - F\{\sigma_2(h)\}.
```

The package finds $`\sigma_1(h)`$ and $`\sigma_2(h)`$ by Newton’s method
in $`m`$, kept inside a bracket by bisection, with the exact derivative
$`\bar S - I`$, and solves this equation for the quantiles of $`H`$.
Near a constant profile, or when $`I`$ is close to 1, $`\bar S`$ and
$`I`$ agree in most of their digits, and their plain difference loses
those digits. To keep them, the package uses the exact identity
$`S_j - I = -(1 - I)\, S_j\{\exp(-\Delta_j) - 1\}`$, evaluated with
[`expm1()`](https://rdrr.io/r/base/Log.html), and subtracts $`I(1 - I)`$
times the mean of the computed $`\Delta_j`$, which rounding moves
slightly away from 0.

### Draws from other software

The diagnostics depend on the draws only through $`\sigma_\tau`$, so
draws of $`\sigma_\tau`$ from any algorithm that targets the same model
and prior give the same posterior summaries, up to Monte Carlo error.
[`reliability_draws()`](https://joonho112.github.io/poolinfo/reference/reliability_draws.md)
computes their Monte Carlo standard errors with the effective sample
sizes of
[`posterior::ess_mean()`](https://mc-stan.org/posterior/reference/ess_mean.html)
([Vehtari et al. 2021](#ref-vehtari2021)). Its site summaries average
over the draws of the mean rather than integrate the mean out. In the
tests, draws from
[`profile_draws()`](https://joonho112.github.io/poolinfo/reference/profile_draws.md)
go to
[`reliability_draws()`](https://joonho112.github.io/poolinfo/reference/reliability_draws.md),
whose posterior means must be within 4.5 Monte Carlo standard errors of
the integrals.

You can compare a brms fit ([Bürkner 2017](#ref-burkner2017)) of the
eight schools with the integrals in the same way. In brms a standard
deviation is bounded below at 0, so a Student t prior of class `sd` is
the half-t of
[`sd_prior()`](https://joonho112.github.io/poolinfo/reference/sd_prior.md),
and an empty prior on the intercept makes the prior on the mean flat.

``` r

library(brms)
fit_brms <- brm(
  estimate | se(se) ~ 1 + (1 | school), data = eight_schools,
  prior = c(prior("", class = Intercept),
            prior(student_t(3, 0, 10), class = sd)),
  chains = 4, iter = 4000, seed = 1, control = list(adapt_delta = 0.99)
)
from_brms <- reliability_draws(fit_brms)
from_brms
exact <- reliability_posterior(estimate, se = se, site = school,
                               data = eight_schools,
                               prior = sd_prior("half_t", scale = 10, df = 3))
exact
```

The posterior means that `from_brms` prints should differ from those of
`exact` by no more than a few of the Monte Carlo standard errors printed
beside them.

### Confidence intervals for the between-site variance

In metafor’s notation $`\sigma_\tau^2`$ is $`\tau^2`$. Suppose an
interval $`C`$ covers $`\sigma_\tau^2`$ with probability $`1 - \alpha`$.
Because $`I`$ increases with $`\sigma_\tau^2`$, the interval between
$`I`$ at the ends of $`C`$ covers $`I`$ with the same probability, and
so does the interval for $`\bar S`$. So does the interval for
$`J_{\mathrm{eff}}`$ when the sampling variances differ; when they are
all equal, $`J_{\mathrm{eff}} = J`$ for every $`\sigma_\tau`$ and its
interval is the single value $`J`$. $`H`$ is not monotone, and
[`map_tau2_interval()`](https://joonho112.github.io/poolinfo/reference/map_tau2_interval.md)
reports its smallest and largest values over $`C`$ instead, using
$`H`$’s peak when the peak lies inside $`C`$. That range contains $`H`$
whenever $`C`$ contains $`\sigma_\tau^2`$, so its coverage is at least
$`1 - \alpha`$.

For a metafor fit, $`C`$ is the Q-profile interval of Viechtbauer
([2007](#ref-viechtbauer2007)) for the model of the fit, moderators
included. With $`Q(\sigma_\tau^2)`$ the weighted residual sum of squares
of the estimates at $`\sigma_\tau^2`$ (the generalized $`Q`$ statistic),
its lower end solves $`Q(\sigma_\tau^2) = \chi^2_{J-p,\,1-\alpha/2}`$
and its upper end $`Q(\sigma_\tau^2) = \chi^2_{J-p,\,\alpha/2}`$, where
$`\chi^2_{J-p,\,q}`$ is the $`q`$ quantile; an end is 0 when $`Q(0)`$ is
already below its quantile.
[`map_tau2_interval()`](https://joonho112.github.io/poolinfo/reference/map_tau2_interval.md)
solves these equations itself, for
$`\log(\sigma_\tau^2 / \mathrm{GM}(v))`$, so that the ends do not depend
on the units of the effects and an end close to 0 has the same relative
precision as one far above it. It does so for every estimator of
metafor, also for `GENQ` and `GENQM` fits, whose own
[`confint()`](https://rdrr.io/r/stats/confint.html) uses another method.
The tests compare the mapped ends with the profiles at the ends of the
interval, and the upper end for $`H`$ with $`H`$’s peak when the peak
lies inside the interval. For the eight schools they also compare the
end found with that of
[`confint()`](https://rdrr.io/r/stats/confint.html) run with a wider
limit. We fit the eight schools with metafor and map its interval.

``` r

reml <- metafor::rma(estimate, sei = se, data = eight_schools,
                     method = "REML")
confint(reml)
#> 
#>        estimate  ci.lb     ci.ub 
#> tau^2    0.0000 0.0000 >100.0000 
#> tau      0.0000 0.0000  >10.0000 
#> I^2(%)   0.0000 0.0000  >42.3703 
#> H^2      1.0000 1.0000   >1.7352
mapped <- map_tau2_interval(reml)
attr(mapped, "interval")
#> estimate    lower    upper 
#>   0.0000   0.0000 270.3672
mapped
#>   functional estimate    lower      upper
#> 1          I 0.000000 0.000000 0.64736618
#> 2       mean 0.000000 0.000000 0.64068042
#> 3      J_eff 6.817468 6.817468 7.76992874
#> 4          H 0.000000 0.000000 0.02836966
```

metafor reports the upper end as above 100, and the interval that is
mapped ends at 270.4. Through $`I`$ the difference is large, since $`I`$
is 0.40 at 100 and 0.65 at 270.4.

## Notation

The paper and the package write some quantities differently. For each
quantity, the table below gives the paper’s symbol, the package’s symbol
or name, and the meaning.

| Paper | Package | Meaning |
|:---|:---|:---|
| $`\hat\tau_j`$, $`v_j`$ | $`\hat\tau_j`$, $`v_j`$, `estimate`, `se2` (`se`$`^2`$) | site estimate and sampling variance |
| $`\sigma`$, $`a = \sigma^2`$ | $`\sigma_\tau`$, $`\sigma_\tau^2`$, `sigma_tau`, `sigma_tau2` | between-site standard deviation and variance |
| $`\tau`$ | $`\mu`$, `effect_mean` | mean of the site effects |
| none | `effect`, a column of [`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md) | posterior mean of a site’s effect |
| none | $`\beta`$, `coef_*` | coefficients of the mean model |
| $`\tilde v`$ | $`\mathrm{GM}(v)`$ | geometric mean of the $`v_j`$ |
| $`g(\sigma)`$ | $`g(\sigma_\tau)`$ | unnormalized posterior density of the between-site standard deviation |
| $`S_j`$ | `reliability`, [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md) | site reliability |
| $`I`$ | `I`, [`kl_centroid()`](https://joonho112.github.io/poolinfo/reference/kl_centroid.md) | informativeness |
| $`\bar S`$ | `mean` | mean reliability |
| $`\mathrm{HM}(S)`$ | `harmonic` | harmonic centroid |
| $`H`$ | `H`, [`kl_dispersion()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md) | dispersion |
| $`V`$, $`\rho_3`$ | `V`, `rho3`, [`logit_moments()`](https://joonho112.github.io/poolinfo/reference/logit_moments.md) | logit variance and third absolute logit moment |
| $`J_{\mathrm{eff}}`$ | `J_eff`, [`concentration_count()`](https://joonho112.github.io/poolinfo/reference/concentration_count.md) | concentration count |
| $`m`$ | `logit_I`, a column of [`profile_path()`](https://joonho112.github.io/poolinfo/reference/profile_path.md) | logit of $`I`$ |
| $`\Delta_j`$ | $`\Delta_j`$ | centered logits, $`\mathrm{logit}\, S_j - m`$ |
| $`R_H`$, $`R_S`$ | `exact` minus `leading`, columns of [`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md) | remainders of Proposition 2 |
| $`\sqrt3\rho_3/108`$, $`\rho_3/48`$ | `bound`, a column of [`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md) | bounds on $`\lvert R_H\rvert`$ and $`\lvert R_S\rvert`$ |
| $`\kappa`$, $`n_j`$ | $`\kappa_j`$, $`n_j`$, `kappa`, `n` | variance constant and site size, $`v_j = \kappa_j/n_j`$ |

## References

Bürkner, Paul-Christian. 2017. “Brms: An R Package for Bayesian
Multilevel Models Using Stan.” *Journal of Statistical Software* 80 (1):
1–28. <https://doi.org/10.18637/jss.v080.i01>.

Harville, David A. 1974. “Bayesian Inference for Variance Components
Using Only Error Contrasts.” *Biometrika* 61 (2): 383–85.
<https://doi.org/10.1093/biomet/61.2.383>.

Lee, JoonHo. 2026. “The Geometry of Informativeness: A Bregman-Centroid
Diagnostic for Partial Pooling in Multisite Trials.” Unpublished
manuscript.

Piessens, Robert, Elise de Doncker-Kapenga, Christoph W. Überhuber, and
David K. Kahaner. 1983. *QUADPACK: A Subroutine Package for Automatic
Integration*. Springer Series in Computational Mathematics. Springer.
<https://doi.org/10.1007/978-3-642-61786-7>.

Vehtari, Aki, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner. 2021. “Rank-Normalization, Folding, and
Localization: An Improved $`\widehat{R}`$ for Assessing Convergence of
MCMC.” *Bayesian Analysis* 16 (2): 667–718.
<https://doi.org/10.1214/20-BA1221>.

Viechtbauer, Wolfgang. 2007. “Confidence Intervals for the Amount of
Heterogeneity in Meta-Analysis.” *Statistics in Medicine* 26 (1): 37–52.
<https://doi.org/10.1002/sim.2514>.

Viechtbauer, Wolfgang. 2010. “Conducting Meta-Analyses in R with the
Metafor Package.” *Journal of Statistical Software* 36 (3): 1–48.
<https://doi.org/10.18637/jss.v036.i03>.
