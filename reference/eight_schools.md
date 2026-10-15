# Eight schools

Estimated effects of coaching programs for the verbal part of the SAT in
eight high schools, with their standard errors (Rubin, 1981). A standard
small example of partial pooling.

## Usage

``` r
eight_schools
```

## Format

A data frame with 8 rows and 3 columns:

- school:

  Label of the school, `"A"` to `"H"`.

- estimate:

  Estimated effect of coaching, in points of the SAT-V.

- se:

  Standard error of the estimate.

## Source

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.), Table 5.2. CRC
Press.

## References

Rubin, D. B. (1981). Estimation in parallel randomized experiments.
*Journal of Educational Statistics, 6*(4), 377-401.
[doi:10.3102/10769986006004377](https://doi.org/10.3102/10769986006004377)

## Examples

``` r
eight_schools
#>   school estimate se
#> 1      A       28 15
#> 2      B        8 10
#> 3      C       -3 16
#> 4      D        7 11
#> 5      E       -1  9
#> 6      F        1 11
#> 7      G       18 10
#> 8      H       12 18
reliability_profile(eight_schools, se = se, sigma_tau = 10)
#> Reliability profile of 8 sites at sigma_tau = 10 (supplied)
#> 
#>   informativeness I            0.404
#>   mean reliability             0.410
#>   dispersion H                 0.0270
#>   logit variance V             0.230
#>   concentration count J_eff    7.5 of 8
#> 
#> Site reliabilities range from 0.24 to 0.55.
```
