# Site allocations under size and budget constraints

Computes allocations under size and budget constraints. The continuous
allocation maximizes the informativeness. With integer sizes and equal
costs the allocation is optimal; with unequal costs the function uses a
greedy rounding rule that need not find the integer optimum.

## Usage

``` r
balanced_allocation(
  J,
  N = NULL,
  cost = NULL,
  budget = NULL,
  lower = NULL,
  upper = NULL,
  integer = TRUE
)
```

## Arguments

- J:

  Number of sites, a positive integer.

- N:

  Total number of individuals, a positive number. Give `N`, for equal
  costs per individual, or `cost` and `budget`, not both. It must be an
  integer when `integer = TRUE` and is then allocated in full.

- cost:

  Costs per individual in each site, positive numbers, one for all sites
  or one per site.

- budget:

  Maximum budget, a positive number, in the units of `cost`.

- lower, upper:

  Optional bounds on the site sizes, numbers, one for all sites or one
  per site.

- integer:

  If `TRUE` (the default), the sizes are integers.

## Value

A numeric vector of `J` site sizes, with the attributes `total` (the
number of individuals) and, with costs, `spent`.

## Details

The informativeness increases with the product of the site sizes when
the constants \\\kappa_j\\ are fixed, so the continuous allocation
maximizes \\\sum_j \log n_j\\ (Lee, 2026, Proposition 4 and Appendix
E.3):

- with equal costs and no bounds, every site gets \\N/J\\.

- with bounds \\L_j \le n_j \le U_j\\, \\n_j = \min\\U_j, \max(L_j,
  t)\\\\ with the level \\t\\ that uses up \\N\\.

- with costs \\c_j\\ and a budget \\B\\, every site gets the same
  spending, \\n_j = B / (J c_j)\\, clipped to the bounds in the same
  way.

With `integer = TRUE`, the continuous solution is rounded down and the
remaining individuals (or budget) go one at a time to the sites where an
individual adds most to \\\sum_j \log n_j\\ per unit of cost. With equal
costs the sizes then differ by at most 1 wherever the bounds allow; with
unequal costs some of the budget can be left over. This greedy rule is
optimal with equal costs, but is a heuristic with unequal costs. For
example, with costs `c(1, 3)` and budget 10 it returns `c(7, 1)`,
whereas `c(4, 2)` has a larger product of site sizes at the same cost.
Bounds are intersected with the positive integers before an integer
allocation is computed; an infeasible request gives an error. A budget
may also be left unused when every site has reached its upper bound.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## See also

[`reliability_design()`](https://joonho112.github.io/poolinfo/reference/reliability_design.md)
to compare a planned allocation with these, and
[`vignette("planning")`](https://joonho112.github.io/poolinfo/articles/planning.md).

## Examples

``` r
balanced_allocation(J = 7, N = 100)
#> [1] 15 15 14 14 14 14 14
#> attr(,"total")
#> [1] 100
balanced_allocation(J = 5, N = 100, upper = c(10, 50, 50, 50, 50))
#> [1] 10 23 23 22 22
#> attr(,"total")
#> [1] 100
balanced_allocation(J = 4, cost = c(1, 1, 2, 4), budget = 400)
#> [1] 100 100  50  25
#> attr(,"total")
#> [1] 275
#> attr(,"spent")
#> [1] 400
```
