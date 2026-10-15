#' Eight schools
#'
#' Estimated effects of coaching programs for the verbal part of the SAT in
#' eight high schools, with their standard errors (Rubin, 1981). A standard
#' small example of partial pooling.
#'
#' @format A data frame with 8 rows and 3 columns:
#' \describe{
#'   \item{school}{Label of the school, `"A"` to `"H"`.}
#'   \item{estimate}{Estimated effect of coaching, in points of the SAT-V.}
#'   \item{se}{Standard error of the estimate.}
#' }
#'
#' @source Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari,
#'   A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.), Table
#'   5.2. CRC Press.
#'
#' @references
#' Rubin, D. B. (1981). Estimation in parallel randomized experiments.
#' *Journal of Educational Statistics, 6*(4), 377-401.
#' \doi{10.3102/10769986006004377}
#'
#' @examples
#' eight_schools
#' reliability_profile(eight_schools, se = se, sigma_tau = 10)
"eight_schools"

#' Example reliability profiles
#'
#' Five sets of sampling variances for 40 sites whose reliability profiles
#' at `sigma_tau = 0.2` all have the informativeness 0.3: equal variances,
#' and four shapes with the logit variance 1.5. The homogeneous, even and
#' concentrated_precise profiles are those of Figure 2 of Lee (2026).
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{profile}{The profile, a factor: `homogeneous`, `even`, `bimodal`,
#'     `concentrated_precise` or `concentrated_noisy`; see [shape_profile()].}
#'   \item{site}{Number of the site, 1 to 40.}
#'   \item{se2, se}{Sampling variance and standard error.}
#'   \item{sigma_tau}{The between-site standard deviation, 0.2.}
#'   \item{reliability}{The site's reliability at that `sigma_tau`.}
#' }
#'
#' @source Built by [shape_profile()], with the script
#'   `data-raw/example_profiles.R` of the source repository.
#'
#' @template ref-lee2026
#' @examples
#' by_profile <- split(example_profiles, example_profiles$profile)
#' sapply(by_profile, function(d) kl_dispersion(d$reliability))
"example_profiles"

#' Online Credit Recovery Study
#'
#' Estimated effects of an online credit recovery course in the 32 sites of
#' the Online Credit Recovery Study (Rickles et al., 2024) analyzed by Lee
#' (2026). In 24 Los Angeles high schools, students who had failed Algebra 1
#' or English 9 in their first year were assigned by lottery to an online
#' course or to the school's usual teacher-directed course, in the summers of
#' 2018 and 2019. A site is a school and a subject: there are 13 Algebra 1
#' sites and 19 English 9 sites, and 8 schools have both.
#'
#' @format A data frame with 32 rows and 18 columns:
#' \describe{
#'   \item{site}{Label of the site, `"S01"` to `"S32"`: the Algebra 1 sites
#'     first, each subject in the order of the study's school identifiers.}
#'   \item{school}{Number of the school, 1 to 24, in the same order. The two
#'     sites of a school have the same number.}
#'   \item{subject}{`"Algebra 1"` or `"English 9"`.}
#'   \item{n, n_treatment, n_control}{Number of students in the site, and of
#'     those assigned to the online course and to the usual course.}
#'   \item{blocks}{Number of randomization blocks in the site.}
#'   \item{estimate, se}{Effect of assignment to the online course on passing
#'     the course, as a difference in proportions, and its standard error:
#'     the coefficient of the assignment in a least squares regression of
#'     passing on the assignment and the randomization blocks, with the HC2
#'     standard error.}
#'   \item{estimate_adjusted, se_adjusted}{The same, with four baseline
#'     covariates added to the regression: the grade 8 state test score in
#'     the subject and the grade 9 GPA, with means imputed for missing
#'     values, and an indicator of a missing value for each.}
#'   \item{estimate_credits, se_credits}{Effect on the cumulative course
#'     credits earned in the subject by the end of the fourth year of high
#'     school, summer terms included, from the same regression as
#'     `estimate`. Credits are the district's course credits; most
#'     semester-long courses earn 5.}
#'   \item{estimate_unadjusted, se_unadjusted}{The difference of the passing
#'     proportions of the two groups, ignoring the blocks, and its usual
#'     standard error.}
#'   \item{log_odds_ratio, se_log_odds_ratio}{Log odds ratio of passing and
#'     its standard error, the square root of the sum of the reciprocals of
#'     the four counts of the site's table of assignment by passing.}
#'   \item{zero_cell}{`TRUE` for a site whose table has a count of 0, in which
#'     case 0.5 is added to each count for the log odds ratio. No site has
#'     one.}
#' }
#'
#' @details
#' The analyses of Lee (2026) are the passing outcome in all 32 sites
#' (`estimate` and `se`) and in each subject, the adjusted passing outcome,
#' the credits, and passing on the log-odds scale. Each is conditional on
#' the standard errors, which are estimated from the students of each site.
#'
#' @source Computed by the package author from the public-use impact analysis
#'   file of the study, `OCRS_publicuse_impactanalysisdata.csv` (American
#'   Institutes for Research, 2024), available from openICPSR, project
#'   202181, \url{https://www.openicpsr.org/openicpsr/project/202181}. The
#'   table contains these site summaries only, not the records of the
#'   students. The public-use file is not part of the package, and its use is
#'   subject to the terms of its release. The script `data-raw/ocrs.R` of the
#'   source repository builds the table from it.
#'
#' @references
#' American Institutes for Research. (2024). *Assessing the efficacy of
#' online credit recovery on student learning and high school graduation:
#' Public-use data file documentation*.
#' \url{https://www.openicpsr.org/openicpsr/project/202181}
#' @template ref-lee2026
#' @references
#' Rickles, J., Clements, M., Brodziak de los Reyes, I., Lachowicz, M., Lin,
#' S., & Heppen, J. (2024). A multisite randomized study of an online
#' learning approach to high school credit recovery: Effects on student
#' experiences and proximal outcomes. *Journal of Research on Educational
#' Effectiveness, 17*(3), 467-490. \doi{10.1080/19345747.2023.2198524}
#' @examples
#' fit <- reliability_posterior(estimate, se = se, site = site, data = ocrs)
#' fit
#'
#' algebra <- subset(ocrs, subject == "Algebra 1")
#' reliability_posterior(estimate, se = se, site = site, data = algebra)
"ocrs"
