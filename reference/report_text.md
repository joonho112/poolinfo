# Sentences for a report

Writes the diagnostics of an analysis as sentences for the results
section of a paper, to be edited into the paper's own prose.

## Usage

``` r
report_text(x, outcome = NULL, digits = 2)
```

## Arguments

- x:

  An object of class `reliability_posterior`, `reliability_draws`,
  `reliability_profile`, `reliability_design` or `profile_comparison`.

- outcome:

  Optional name of the outcome, a string, as in `"course passing"`.

- digits:

  Number of decimals of the informativeness, the mean reliability and
  the probability, a positive integer; the dispersion has one more.

## Value

A character vector with one sentence per element; for a comparison, the
sentences of each analysis in turn.

## Details

The sentences open with the scale of the effects and the estimand at
each site when the fit records them (`scale_label` and `target_label`).
For a posterior fit or draws, the next sentence states the working model
and the prior; for a profile, the between-site standard deviation and
where it came from; and for a design, the planned allocation and its
comparison with balance or equal spending. When the threshold is 1/2, a
sentence says what that value of the informativeness means: the
between-site variance equals the geometric mean of the sampling
variances. The last sentence states that the summaries are conditional
on the supplied standard errors. The sentences give numbers and
intervals without describing them as low or high.

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, data = eight_schools,
                             scale_label = "SAT-V point",
                             target_label = "the effect of coaching on the SAT-V")
report_text(fit, outcome = "the SAT-V")
#> [1] "The effects are on the SAT-V point scale."                                                                                                                                                                                                                         
#> [2] "The estimand at each site is the effect of coaching on the SAT-V."                                                                                                                                                                                                 
#> [3] "Under the normal working model with flat priors on the mean and the between-site standard deviation, the informativeness of the 8 sites for the SAT-V was 0.23 (95% posterior interval 0.00 to 0.75), and the posterior probability that it exceeded 0.5 was 0.14."
#> [4] "At an informativeness of 0.5, the between-site variance equals the geometric mean of the sampling variances."                                                                                                                                                      
#> [5] "The mean reliability was 0.23 (0.00 to 0.74)."                                                                                                                                                                                                                     
#> [6] "The dispersion of the site reliabilities was H = 0.014 (0.000 to 0.028), with logit variance V = 0.23, and the concentration count was 7.2 (at most 8, the number of sites)."                                                                                      
#> [7] "These summaries are conditional on the supplied standard errors, and each site counts equally in them, whatever its size."                                                                                                                                         
```
