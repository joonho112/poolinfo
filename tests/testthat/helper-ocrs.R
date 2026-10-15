# The six OCRS analyses of Lee (2026): the inputs from the data set ocrs and
# reference posterior summaries computed for the paper, kept in
# cases/ocrs-reference.rds.
ocrs_analyses <- function() {
  reference <- readRDS(testthat::test_path("cases", "ocrs-reference.rds"))$analyses
  algebra <- ocrs$subject == "Algebra 1"
  inputs <- list(
    pooled = list(y = ocrs$estimate, v = ocrs$se^2),
    algebra = list(y = ocrs$estimate[algebra], v = ocrs$se[algebra]^2),
    english = list(y = ocrs$estimate[!algebra], v = ocrs$se[!algebra]^2),
    pooled_adj = list(y = ocrs$estimate_adjusted, v = ocrs$se_adjusted^2),
    credits = list(y = ocrs$estimate_credits, v = ocrs$se_credits^2),
    logodds = list(y = ocrs$log_odds_ratio, v = ocrs$se_log_odds_ratio^2)
  )
  stopifnot(identical(names(inputs), names(reference)))
  Map(c, inputs, reference)
}

# The profile of the passing outcome at sigma_tau = 0.188, the posterior mean
# of sigma_tau for that outcome, with its reference values.
ocrs_profile_reference <- function() {
  readRDS(testthat::test_path("cases", "ocrs-reference.rds"))$profile
}

# Six simulated data sets with their reference posterior summaries.
read_boundary_cases <- function() readRDS(testthat::test_path("cases", "boundary-cases.rds"))
