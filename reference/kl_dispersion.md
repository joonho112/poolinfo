# Dispersion of a reliability profile

`kl_dispersion()` computes the dispersion \\H\\ of the site
reliabilities around their informativeness \\I\\, and
`kl_contributions()` the term of each site, whose mean is \\H\\.

## Usage

``` r
kl_dispersion(S)

kl_contributions(S)
```

## Arguments

- S:

  Site reliabilities, a numeric vector of numbers strictly between 0 and
  1, as returned by
  [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md).

## Value

`kl_dispersion()` returns \\H\\, a number of at least 0.
`kl_contributions()` returns a numeric vector with one term for each
site, with the names of `S`.

## Details

The dispersion is the minimum of the average divergence that defines
\\I\\, \$\$H = \frac{1}{J}\sum_j \mathrm{KL}\\\mathrm{Ber}(I) \\
\mathrm{Ber}(S_j)\\ = \frac{1}{J}\sum_j A(m + \Delta_j) - A(m),\$\$
where \\A(t) = \log(1 + e^t)\\, \\m\\ is the logit of \\I\\ and
\\\Delta_j\\ the centered logit of site \\j\\ (Lee, 2026, Theorem 2 and
Proposition 2). It is 0 only when all the reliabilities are equal.

A site's term is the divergence of \\I\\ from the site's reliability.
Removing the site would also move \\I\\, so the term is not the change
in \\H\\ without that site.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`logit_moments()`](https://joonho112.github.io/poolinfo/reference/logit_moments.md)
and
[`local_terms()`](https://joonho112.github.io/poolinfo/reference/local_terms.md)
for the logit variance \\V\\ and the second-order approximation \\H
\approx I(1-I)V/2\\.

## Examples

``` r
S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
kl_dispersion(S)
#> [1] 0.07569242
kl_contributions(S)
#> [1] 0.134361601 0.009963631 0.010239076 0.148205370
mean(kl_contributions(S))
#> [1] 0.07569242
```
