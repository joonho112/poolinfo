# Reliability profile of a planned allocation

Computes the reliability profile of a planned allocation of individuals
to sites and compares it with the balanced allocation of the same total,
or with equal spending when the costs per individual differ.

## Usage

``` r
reliability_design(
  n,
  sigma_tau,
  kappa = NULL,
  p = 0.5,
  R2 = 0,
  var_outcome = 1,
  cost = NULL
)
```

## Arguments

- n:

  Planned numbers of individuals in the sites, a numeric vector of
  positive numbers.

- sigma_tau:

  Planned between-site standard deviation of the effects, a positive
  number, or a numeric vector of them for a range of scenarios. The
  printed comparison uses the first.

- kappa:

  Optional constants \\\kappa_j\\ with \\v_j = \kappa_j / n_j\\, one for
  all sites or one per site. By default they come from `p`, `R2` and
  `var_outcome`, as in
  [`neyman_variance()`](https://joonho112.github.io/poolinfo/reference/neyman_variance.md).

- p:

  Proportion of each site's individuals assigned to treatment, numbers
  strictly between 0 and 1, one for all sites or one per site.

- R2:

  Proportion of the outcome variance explained by covariates, numbers in
  \\\[0, 1)\\, one for all sites or one per site.

- var_outcome:

  Variance of the outcome within each arm of a site (\\s_j^2\\ in
  equation 2 of Lee, 2026), positive numbers, one for all sites or one
  per site.

- cost:

  Optional costs per individual in each site, positive numbers. The
  comparison is then with equal spending of the same budget.

## Value

An object of class `reliability_design`, a list whose main components
are

- `n`, `kappa`, `se2`, `sigma_tau`, `cost`: the inputs, with the
  sampling variances `se2`.

- `planned`, `balanced`: `reliability_profile` objects of the planned
  allocation and of the balanced one (or equal spending) at the first
  `sigma_tau`.

- `balanced_n`: the sizes of the balanced allocation or of equal
  spending.

- `efficiency`: the ratio of the geometric to the arithmetic mean of the
  planned sizes.

- `local`: the approximation of Appendix E.2, a list with the
  coefficient of variation of the sizes `cv`, the logarithm of the
  efficiency `exact`, its approximation `leading` (\\-\mathrm{CV}^2/2\\)
  and the bound on the remainder `bound` (`NA` when a size is twice the
  mean or more).

- `scenarios`: a data frame with the informativeness of both allocations
  at each `sigma_tau`.

[`print()`](https://rdrr.io/r/base/print.html),
[`profile_summary()`](https://joonho112.github.io/poolinfo/reference/profile_summary.md)
(with a column `allocation`),
[`site_table()`](https://joonho112.github.io/poolinfo/reference/site_table.md)
(with the columns `n`, `kappa` and `n_balanced`),
[`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
(the planned profile),
[`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md)
and
[`plot_design()`](https://joonho112.github.io/poolinfo/reference/plot_design.md)
take it.

## Details

With a fixed number of sites \\J\\, a fixed total \\N\\ and a common
\\\kappa\\, the balanced allocation \\n_j = N/J\\ maximizes the
informativeness and the mean reliability, and minimizes the benchmark
risk; there \\I = \sigma\_\tau^2 / (\sigma\_\tau^2 + \kappa J / N)\\
(Lee, 2026, Proposition 4). The informativeness depends on the sizes
through their geometric mean \\\tilde n\\: with \\u_j = n_j / \bar n -
1\\ and \\\mathrm{CV}\_J^2 = J^{-1}\sum_j u_j^2\\, \\\log(\tilde n /
\bar n) = -\mathrm{CV}\_J^2/2 + R\\ with \\\|R\| \le J^{-1}\sum_j
\|u_j\|^3 / \\3(1 - r)^3\\\\ when every \\\|u_j\| \le r \< 1\\ (Appendix
E.2). The component `local` contains this approximation and its bound.

With unequal costs \\c_j\\, the allocation that maximizes the
informativeness for a budget \\B\\ spends the same on every site, \\n_j
= B / (J c_j)\\ (Appendix E.3). With constants \\\kappa_j\\ that differ,
balance still maximizes the informativeness, because \\\mathrm{GM}(v) =
\mathrm{GM}(\kappa)/\mathrm{GM}(n)\\, but it need not make the
reliabilities equal or maximize the mean reliability (Appendix E.3).

These results concern the site profile at known inputs: the
informativeness, the mean reliability and the benchmark risk, not ranks
or the distribution of the effects.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`balanced_allocation()`](https://joonho112.github.io/poolinfo/reference/balanced_allocation.md)
and
[`solve_design()`](https://joonho112.github.io/poolinfo/reference/solve_design.md)
for allocations and sizes, and
[`vignette("planning")`](https://joonho112.github.io/poolinfo/articles/planning.md)
for a worked plan.

## Examples

``` r
n <- round(exp(qnorm(seq(0.5 / 40, 1 - 0.5 / 40, length.out = 40), 3.8, 0.47)))
design <- reliability_design(n, sigma_tau = 0.2, kappa = 4)
design
#> Planned allocation: 40 sites, 1,990 individuals, kappa = 4, sigma_tau = 0.2
#> Site sizes from 16 to 128 (CV 0.48)
#> 
#>                              planned    balanced
#>   informativeness I          0.309      0.332
#>   mean reliability           0.317      0.332
#>   dispersion H               0.0225     0.0000
#>   concentration count J_eff  36.6       40.0
#>   geometric mean size        44.7       49.8
#> 
#> Balance maximizes I and the mean reliability at fixed J and N when every site
#> has the same kappa and the same cost per individual (Lee, 2026).

# The informativeness of both allocations over a range of sigma_tau
reliability_design(n, sigma_tau = c(0.1, 0.2, 0.3), kappa = 4)$scenarios
#>   sigma_tau I_planned I_balanced
#> 1       0.1 0.1006091  0.1106170
#> 2       0.2 0.3091319  0.3322204
#> 3       0.3 0.5016873  0.5281628
```
