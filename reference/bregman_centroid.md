# Bregman centroids and dispersions of a reliability profile

`bregman_centroid()` computes the centroid of the site reliabilities
under one of three Bregman divergences, on either side of the
divergence, and `bregman_dispersion()` the average divergence at that
centroid.

## Usage

``` r
bregman_centroid(
  S,
  potential = c("kl", "squared", "burg"),
  side = c("first", "second")
)

bregman_dispersion(
  S,
  potential = c("kl", "squared", "burg"),
  side = c("first", "second")
)
```

## Arguments

- S:

  Site reliabilities, a numeric vector of numbers strictly between 0 and
  1, as returned by
  [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md).

- potential:

  The convex function that defines the divergence: one of `"kl"` (the
  default), `"squared"` or `"burg"`. See Details.

- side:

  Which argument of the divergence the centroid takes: `"first"` (the
  default) or `"second"`.

## Value

A number: the centroid, between 0 and 1, or the dispersion, of at least
0.

## Details

A strictly convex function \\\phi\\ defines the divergence \\D\_\phi(x,
y) = \phi(x) - \phi(y) - \phi'(y)(x - y)\\. The first-argument centroid
minimizes \\\sum_j D\_\phi(c, S_j)\\ and is \\(\phi')^{-1}\\J^{-1}\sum_j
\phi'(S_j)\\\\; the second-argument centroid minimizes \\\sum_j
D\_\phi(S_j, c)\\ and is the mean reliability for every \\\phi\\ (Lee,
2026, Table 2 and Appendix B):

|             |                              |                         |
|-------------|------------------------------|-------------------------|
| `potential` | \\\phi(x)\\                  | first-argument centroid |
| `"squared"` | \\x^2/2\\                    | mean reliability        |
| `"kl"`      | \\x\log x + (1-x)\log(1-x)\\ | informativeness \\I\\   |
| `"burg"`    | \\-\log x\\                  | harmonic mean           |

The dispersion is the minimum average divergence. On the first side of
`"kl"` it is the dispersion \\H\\ of
[`kl_dispersion()`](https://joonho112.github.io/poolinfo/reference/kl_dispersion.md);
on the second side it is the mean Bernoulli divergence of the
reliabilities from their mean. The two need not agree: for the
reliabilities 0.2 and 0.8 both centroids are 0.5, but the dispersions
are \\\log(5/4) \approx 0.2231\\ and about 0.1927. For `"burg"` the
dispersions are the logarithm of the ratio of the geometric to the
harmonic mean (first side) and of the arithmetic to the geometric mean
(second side); for `"squared"` both are half the variance of the
reliabilities, with divisor \\J\\.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
S <- c(0.2, 0.8)
bregman_centroid(S, "kl")
#> [1] 0.5
bregman_dispersion(S, "kl", side = "first")
#> [1] 0.2231436
bregman_dispersion(S, "kl", side = "second")
#> [1] 0.1927448

S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
sapply(c("squared", "kl", "burg"), function(p) bregman_centroid(S, p))
#>   squared        kl      burg 
#> 0.5619231 0.5714286 0.4961240 
```
