# Table of the diagnostics for a report

Formats the diagnostics of one analysis or of a comparison as a table in
the layout of Table 5 of Lee (2026): the informativeness, the mean
reliability and the dispersion with their posterior intervals, the
concentration count and the posterior probability that the
informativeness exceeds the threshold.

## Usage

``` r
report_table(x, digits = 2, format = c("data.frame", "markdown", "latex"))
```

## Arguments

- x:

  An object of class `reliability_posterior`, `reliability_draws`,
  `reliability_profile` or `profile_comparison`.

- digits:

  Number of decimals of the informativeness, the mean reliability and
  the probability, a positive integer; the dispersion has one more and
  the concentration count 1.

- format:

  `"data.frame"` (the default) for a data frame of numbers, `"markdown"`
  or `"latex"` for text in the style of the paper.

## Value

For `format = "data.frame"`, a data frame with one row per analysis and
the columns `analysis`, `J`, `I`, `I_lower`, `I_upper`, `mean`, `H`,
`H_lower`, `H_upper`, `J_eff` and `prob_I_above`, and the note in its
attribute `note`. For the text formats, a character string.

## Details

In the text formats the leading zeros of the numbers between 0 and 1 are
dropped and intervals are in brackets, as in the paper. A table of one
analysis has no column for its name; the names of a comparison label its
rows. The LaTeX table uses the rules of the booktabs package. Every
format ends with a note that states the working model, the prior, the
scale and the estimand when the fits record them, the conditioning on
the standard errors, and that each site counts equally in the
diagnostics.

## References

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools)
report_table(fit)
#>   analysis J         I      I_lower   I_upper      mean          H      H_lower
#> 1 analysis 8 0.2252945 0.0003997828 0.7451728 0.2297284 0.01428773 4.407725e-05
#>      H_upper    J_eff prob_I_above
#> 1 0.02832835 7.170045    0.1379141
cat(report_table(fit, format = "markdown"))
#> | J | I [95% interval] | Mean reliability | H [95% interval] | J_eff | Pr(I > 0.5) |
#> |:--|:--|:--|:--|:--|:--|
#> | 8 | .23 [.00, .75] | .23 | .014 [.000, .028] | 7.2 | .14 |
#> 
#> *Note.* Each estimate is modeled as normal around its site's effect with a known sampling variance, and the site effects as independent and normal across sites. Flat priors on the mean and the between-site standard deviation. Summaries are conditional on the supplied standard errors, and each site counts equally in the diagnostics, whatever its size.
```
