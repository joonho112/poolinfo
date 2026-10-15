# Online Credit Recovery Study

Estimated effects of an online credit recovery course in the 32 sites of
the Online Credit Recovery Study (Rickles et al., 2024) analyzed by Lee
(2026). In 24 Los Angeles high schools, students who had failed Algebra
1 or English 9 in their first year were assigned by lottery to an online
course or to the school's usual teacher-directed course, in the summers
of 2018 and 2019. A site is a school and a subject: there are 13 Algebra
1 sites and 19 English 9 sites, and 8 schools have both.

## Usage

``` r
ocrs
```

## Format

A data frame with 32 rows and 18 columns:

- site:

  Label of the site, `"S01"` to `"S32"`: the Algebra 1 sites first, each
  subject in the order of the study's school identifiers.

- school:

  Number of the school, 1 to 24, in the same order. The two sites of a
  school have the same number.

- subject:

  `"Algebra 1"` or `"English 9"`.

- n, n_treatment, n_control:

  Number of students in the site, and of those assigned to the online
  course and to the usual course.

- blocks:

  Number of randomization blocks in the site.

- estimate, se:

  Effect of assignment to the online course on passing the course, as a
  difference in proportions, and its standard error: the coefficient of
  the assignment in a least squares regression of passing on the
  assignment and the randomization blocks, with the HC2 standard error.

- estimate_adjusted, se_adjusted:

  The same, with four baseline covariates added to the regression: the
  grade 8 state test score in the subject and the grade 9 GPA, with
  means imputed for missing values, and an indicator of a missing value
  for each.

- estimate_credits, se_credits:

  Effect on the cumulative course credits earned in the subject by the
  end of the fourth year of high school, summer terms included, from the
  same regression as `estimate`. Credits are the district's course
  credits; most semester-long courses earn 5.

- estimate_unadjusted, se_unadjusted:

  The difference of the passing proportions of the two groups, ignoring
  the blocks, and its usual standard error.

- log_odds_ratio, se_log_odds_ratio:

  Log odds ratio of passing and its standard error, the square root of
  the sum of the reciprocals of the four counts of the site's table of
  assignment by passing.

- zero_cell:

  `TRUE` for a site whose table has a count of 0, in which case 0.5 is
  added to each count for the log odds ratio. No site has one.

## Source

Computed by the package author from the public-use impact analysis file
of the study, `OCRS_publicuse_impactanalysisdata.csv` (American
Institutes for Research, 2024), available from openICPSR, project
202181, <https://www.openicpsr.org/openicpsr/project/202181>. The table
contains these site summaries only, not the records of the students. The
public-use file is not part of the package, and its use is subject to
the terms of its release. The script `data-raw/ocrs.R` of the source
repository builds the table from it.

## Details

The analyses of Lee (2026) are the passing outcome in all 32 sites
(`estimate` and `se`) and in each subject, the adjusted passing outcome,
the credits, and passing on the log-odds scale. Each is conditional on
the standard errors, which are estimated from the students of each site.

## References

American Institutes for Research. (2024). *Assessing the efficacy of
online credit recovery on student learning and high school graduation:
Public-use data file documentation*.
<https://www.openicpsr.org/openicpsr/project/202181>

Lee, J. (2026). The geometry of informativeness: A Bregman-centroid
diagnostic for partial pooling in multisite trials. Unpublished
manuscript.

Rickles, J., Clements, M., Brodziak de los Reyes, I., Lachowicz, M.,
Lin, S., & Heppen, J. (2024). A multisite randomized study of an online
learning approach to high school credit recovery: Effects on student
experiences and proximal outcomes. *Journal of Research on Educational
Effectiveness, 17*(3), 467-490.
[doi:10.1080/19345747.2023.2198524](https://doi.org/10.1080/19345747.2023.2198524)

## Examples

``` r
fit <- reliability_posterior(estimate, se = se, site = site, data = ocrs)
fit
#> Reliability profile of 32 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#>                              mean     95% interval
#>   informativeness I          0.671    0.497    0.812
#>   mean reliability           0.663    0.498    0.801
#>   dispersion H               0.0253   0.0185   0.0289
#>   concentration count J_eff  31.2     30.3     31.7     of 32
#>   logit variance V           0.237    (fixed by the standard errors)
#>   sigma_tau                  0.188    0.127    0.266
#>   effect mean                -0.130   -0.211   -0.0498
#> 
#>   Pr(I > 0.5)                0.973
#>   I at E(sigma_tau^2)        0.691

algebra <- subset(ocrs, subject == "Algebra 1")
reliability_posterior(estimate, se = se, site = site, data = algebra)
#> Reliability profile of 13 sites (posterior)
#> Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
#>        sites independent
#> Prior: flat on the mean and on sigma_tau
#> Summaries are conditional on the supplied standard errors.
#> Each site counts equally in the diagnostics, whatever its size.
#> 
#>                              mean     95% interval
#>   informativeness I          0.629    0.294    0.865
#>   mean reliability           0.625    0.300    0.859
#>   dispersion H               0.0164   0.0086   0.0196
#>   concentration count J_eff  12.7     12.0     13.0     of 13
#>   logit variance V           0.158    (fixed by the standard errors)
#>   sigma_tau                  0.189    0.0862   0.339
#>   effect mean                -0.0806  -0.213   0.0529
#> 
#>   Pr(I > 0.5)                0.818
#>   I at E(sigma_tau^2)        0.691
```
