test_that("single numbers are checked against their range", {
  expect_identical(.rp_check_level(0.95), 0.95)
  expect_null(.rp_check_level(NULL, null_ok = TRUE))
  for (bad in list(0, 1, -0.1, NA_real_, Inf, "0.9", c(0.5, 0.9), NULL)) {
    expect_rp_error(.rp_check_level(bad), regexp = "between 0 and 1")
  }
  expect_identical(.rp_check_sigma_tau(0), 0)
  expect_identical(.rp_check_sigma_tau(0.2), 0.2)
  expect_rp_error(.rp_check_sigma_tau(-0.1), regexp = "at least 0")
  expect_rp_error(.rp_abort_missing_sigma_tau(), regexp = "missing")
  expect_identical(.rp_check_whole(3, "n"), 3L)
  expect_rp_error(.rp_check_whole(2.5, "n"), regexp = "an integer")
  expect_rp_error(.rp_check_whole(0, "n"), regexp = "at least 1")
  expect_rp_error(.rp_check_whole(9, "n", upper = 5), regexp = "between 1 and 5")
  expect_true(.rp_check_flag(TRUE, "x"))
  expect_rp_error(.rp_check_flag(NA, "x"), regexp = "TRUE or FALSE")
  expect_rp_error(.rp_check_string(1, "label"), regexp = "single string")
  expect_null(.rp_check_string(NULL, "label"))
})

test_that("a choice is one of its options, and the whole set means the first", {
  choices <- c("even", "bimodal")
  expect_identical(.rp_match(choices, choices, "shape"), "even")
  expect_identical(.rp_match("bimodal", choices, "shape"), "bimodal")
  err <- expect_rp_error(.rp_match("bi", choices, "shape"), regexp = "one of \"even\", \"bimodal\"")
  expect_match(conditionMessage(err), "You passed \"bi\"", fixed = TRUE)
})

test_that("unused arguments in the dots are named", {
  expect_true(.rp_check_dots_empty())
  err <- expect_rp_error(.rp_check_dots_empty(levle = 0.9))
  expect_match(conditionMessage(err), "`levle`", fixed = TRUE)
  err <- expect_rp_error(.rp_check_dots_empty(3, 4))
  expect_match(conditionMessage(err), "2 unnamed arguments", fixed = TRUE)
})

test_that("vectors must hold finite numbers, and positive ones where asked", {
  expect_identical(.rp_check_vector(1:3, "x"), c(1, 2, 3))
  expect_identical(.rp_check_vector(matrix(1:3, ncol = 1), "x"), c(1, 2, 3))
  expect_rp_error(.rp_check_vector(matrix(1:4, 2), "x"), regexp = "numeric vector")
  expect_rp_error(.rp_check_vector("a", "x"), regexp = "numeric vector")
  expect_rp_error(.rp_check_vector(numeric(), "x"), regexp = "at least 1 number")
  err <- expect_rp_error(.rp_check_vector(c(1, NA, 3, Inf), "x"), regexp = "finite")
  expect_match(conditionMessage(err), "Positions 2 and 4 are", fixed = TRUE)
  err <- expect_rp_error(.rp_check_vector(c(1, 0, 3), "se2", positive = TRUE), regexp = "positive")
  expect_match(conditionMessage(err), "Position 2 is 0.", fixed = TRUE)
  expect_match(conditionMessage(err), "reliability 1", fixed = TRUE)
  err <- expect_rp_error(.rp_check_vector(c(-1, 2), "se", positive = TRUE))
  expect_match(conditionMessage(err), "Position 1 is negative.", fixed = TRUE)
})

test_that("exactly one of se and se2 is given", {
  expect_equal(.rp_variances(c(a = 0.1, b = 0.2), NULL), c(a = 0.01, b = 0.04), tolerance = 1e-15)
  expect_identical(.rp_variances(NULL, c(0.01, 0.04)), c(0.01, 0.04))
  expect_rp_error(.rp_variances(0.1, 0.01), regexp = "both supplied")
  expect_rp_error(.rp_variances(NULL, NULL), regexp = "standard errors are missing")
  expect_rp_error(.rp_variances(c(0.1, 0), NULL), regexp = "`se` must hold positive")
  expect_rp_error(.rp_variances(NULL, 0.01, min_length = 2L), regexp = "at least 2 numbers")
})

test_that("reliabilities lie strictly between 0 and 1", {
  expect_identical(.rp_check_reliabilities(c(x = 0.2, y = 0.8)), c(x = 0.2, y = 0.8))
  err <- expect_rp_error(.rp_check_reliabilities(c(0.2, 1)), regexp = "strictly between 0 and 1")
  expect_match(conditionMessage(err), "no finite logit", fixed = TRUE)
  expect_rp_error(.rp_check_reliabilities(c(0.2, 1.5)))
  expect_rp_error(.rp_check_reliabilities(c(0.2, NA)), regexp = "finite")
})

test_that("site labels are distinct and one for each site", {
  expect_null(.rp_check_site(NULL, 3))
  expect_identical(.rp_check_site(factor(c("a", "b")), 2), c("a", "b"))
  expect_rp_error(.rp_check_site(c("a", "b"), 3), regexp = "one label for each of the 3 sites")
  expect_rp_error(.rp_check_site(c("a", "a"), 2), regexp = "distinct")
  expect_rp_error(.rp_check_site(c("a", NA), 2), regexp = "missing")
  expect_identical(.rp_site_labels(NULL, 3), c("1", "2", "3"))
  expect_identical(.rp_site_labels(NULL, 2, names = c("x", "y")), c("x", "y"))
  expect_identical(.rp_site_labels(NULL, 2, names = c("x", "x")), c("1", "2"))
})

test_that("an unnamed second argument is refused rather than taken as standard errors", {
  y <- eight_schools$estimate
  v <- eight_schools$se^2
  expect_rp_error(reliability_posterior(y, v), regexp = "ambiguous")
  expect_rp_error(reliability_posterior(estimate, se, data = eight_schools), regexp = "ambiguous")
  expect_rp_error(reliability_draws(abs(stats::rnorm(50, 8, 4)), v), regexp = "ambiguous")
  expect_rp_error(site_reliability(v, sigma_tau = 10), regexp = "ambiguous")
  expect_rp_error(typical_variance(v), regexp = "ambiguous")
  expect_rp_error(benchmark_risk(v, sigma_tau = 10), regexp = "ambiguous")
  expect_rp_error(reliability_posterior(estimate, se = se, sites = school, data = eight_schools),
                  regexp = "no further arguments")
})

test_that("profiles and draws take their columns from data", {
  a <- reliability_profile(se = se, site = school, estimate = estimate, data = eight_schools,
                           sigma_tau = 10)
  b <- reliability_profile(eight_schools, se = se, site = school, estimate = estimate, sigma_tau = 10)
  expect_equal(a$S, b$S)
  expect_equal(site_table(a), site_table(b))
  withr::local_seed(94)
  s <- abs(stats::rnorm(300, 8, 4))
  c1 <- reliability_draws(s, se = se, site = school, data = eight_schools)
  c2 <- reliability_draws(s, se = eight_schools$se, site = eight_schools$school)
  expect_equal(c1$summary, c2$summary)
  expect_rp_error(reliability_draws(s, se = se, data = list(se = 1)), regexp = "data frame")
})
