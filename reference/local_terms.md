# Second-order terms of the dispersion and of the centroid gap

Compares the dispersion \\H\\ and the gap between the mean reliability
and the informativeness with their second-order approximations, and
bounds the remainders.

## Usage

``` r
local_terms(S)
```

## Arguments

- S:

  Site reliabilities, a numeric vector of numbers strictly between 0 and
  1, as returned by
  [`site_reliability()`](https://joonho112.github.io/poolinfo/reference/site_reliability.md).

## Value

A data frame with two rows, `H` and `gap` (the mean reliability minus
\\I\\), and the columns

- `quantity`: `"H"` or `"gap"`.

- `exact`: the value of the quantity.

- `leading`: its second-order term.

- `bound`: the bound on the absolute remainder.

- `holds`: `TRUE` if the absolute leading term is larger than the bound,
  `NA` if both are 0, as when all the reliabilities are equal.

## Details

Proposition 2 of Lee (2026) gives, with the logit variance \\V\\ and the
third absolute logit moment \\\rho_3\\ of
[`logit_moments()`](https://joonho112.github.io/poolinfo/reference/logit_moments.md),
\$\$H = \tfrac12 I(1-I)V + R_H, \qquad \|R_H\| \le
\frac{\sqrt3}{108}\rho_3,\$\$ \$\$\bar S - I = \tfrac12 I(1-I)(1-2I)V +
R_S, \qquad \|R_S\| \le \frac{\rho_3}{48}.\$\$ When the leading term is
larger than the bound in absolute value, the exact value has the sign of
the leading term. For the gap, that sign is the sign of \\1 - 2I\\, so
the mean reliability is above \\I\\ when \\I \< 1/2\\. When the bound is
the larger, the remainder can decide the sign, and the gap can have the
other sign; the `exact` column always gives the value itself.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
local_terms(S)
#>   quantity        exact     leading      bound holds
#> 1        H  0.075692420  0.07896183 0.01082354  TRUE
#> 2      gap -0.009505495 -0.01128026 0.01406019 FALSE

# A skewed profile of 25 sites with I = 0.3 and V = 2 whose mean
# reliability is below I, against the sign of the leading term.
d <- c(sqrt(48), rep(-sqrt(1 / 12), 24))
local_terms(plogis(qlogis(0.3) + d))
#>   quantity       exact leading     bound holds
#> 1        H  0.15398567   0.210 0.2137037 FALSE
#> 2      gap -0.02675316   0.084 0.2776093 FALSE
```
