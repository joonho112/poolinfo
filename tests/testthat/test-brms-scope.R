# These stand-ins carry the parts of a brms fit that the package reads;
# they are not sampled fits.
brms_metadata <- function(formula = y | se(s) ~ 1 + (1 | site),
                          family = stats::gaussian()) {
  structure(list(
    formula = list(formula = formula), family = family,
    data = data.frame(y = c(-0.3, 0.1, 0.2, 0.4), s = rep(0.1, 4),
                      x = c(-1, 0, 1, 2), site = letters[1:4], other = LETTERS[1:4])
  ), class = "brmsfit")
}

test_that("a supported brms model is recognized", {
  x <- brms_metadata()
  parts <- .rp_brms_parts(x)
  expect_identical(parts$response, "y")
  expect_identical(parts$se, "s")
  expect_identical(parts$group, "site")
  expect_identical(parts$mean_variable, "b_Intercept")
  expect_identical(.rp_brms_parts(brms_metadata(y | se(s, sigma = FALSE) ~ 1 + (1 | site))), parts)
  expect_identical(.rp_brms_parts(brms_metadata(y | se(s, FALSE) ~ 1 + (1 || site))), parts)
  expect_identical(.rp_brms_parts(brms_metadata(y | se(s) ~ (1 | site)))$mean_variable, "b_Intercept")
})

test_that("brms models outside the scope are refused before any draws are read", {
  testthat::local_mocked_bindings(.rp_extract_draws = function(...) stop("draw conversion was reached"))
  unsupported <- list(
    brms_metadata(y | se(s, sigma = TRUE) ~ 1 + (1 | site)),
    brms_metadata(family = list(family = "student", link = "identity")),
    brms_metadata(family = stats::gaussian(link = "log")),
    brms_metadata(y | se(s) ~ 1 + (1 | site) + (1 | other)),
    brms_metadata(y | se(s) ~ 1 + (1 + x | site)),
    brms_metadata(y | se(s) ~ 1 + (1 | gr(site))),
    brms_metadata(log(y) | se(s) ~ 1 + (1 | site)),
    brms_metadata(y | se(sqrt(s)) ~ 1 + (1 | site)),
    brms_metadata(y | se(s) + weights(x) ~ 1 + (1 | site)),
    brms_metadata(y | se(s) ~ 1 + (1 | site) + ar(time = x)),
    brms_metadata(y | se(s) ~ 1 + (1 | site) + brms::ar(time = x)),
    brms_metadata(y | se(s) ~ 1 + (1 | site) + gp(x)),
    brms_metadata(y | se(s) ~ 1 + (1 | site) + s(x)),
    brms_metadata(y ~ 1 + (1 | site))
  )
  for (x in unsupported) {
    expect_rp_error(reliability_draws(x), class = "rp_scope_error")
  }
  duplicate <- brms_metadata()
  duplicate$data$site[4] <- duplicate$data$site[1]
  expect_rp_error(reliability_draws(duplicate), class = "rp_scope_error", regexp = "one observation")
  distributional <- brms_metadata()
  distributional$formula$pforms <- list(sigma = sigma ~ x)
  expect_rp_error(reliability_draws(distributional), class = "rp_scope_error")
  autocor <- brms_metadata()
  attr(autocor$formula$formula, "autocor") <- ~ ar(time)
  expect_rp_error(reliability_draws(autocor), class = "rp_scope_error")
})

test_that("a supported brms model gives draws and site summaries", {
  testthat::local_mocked_bindings(.rp_extract_draws = function(x, variable, mean) {
    expect_identical(variable, "sd_site__Intercept")
    expect_identical(mean, "b_Intercept")
    list(sigma = c(0.15, 0.2, 0.25, 0.3), mean = c(-0.1, 0, 0.1, 0.2),
         chains = rep(1L, 4), has_chains = TRUE)
  })
  x <- reliability_draws(brms_metadata())
  expect_equal(x$se2, rep(0.01, 4))
  expect_equal(x$estimate, c(-0.3, 0.1, 0.2, 0.4))
  expect_true(all(c("effect", "effect_sd") %in% names(site_table(x))))
})

test_that("fixed-moderator brms models supply variance-based diagnostics only", {
  testthat::local_mocked_bindings(.rp_extract_draws = function(x, variable, mean) {
    expect_null(mean)
    list(sigma = c(0.15, 0.2, 0.25, 0.3), mean = NULL,
         chains = rep(1L, 4), has_chains = TRUE)
  })
  metadata <- brms_metadata(y | se(s) ~ 1 + x + (1 | site))
  x <- reliability_draws(metadata)
  direct <- reliability_draws(c(0.15, 0.2, 0.25, 0.3), se = rep(0.1, 4))
  expect_equal(x$summary, direct$summary)
  expect_null(x$estimate)
  expect_false("effect" %in% names(site_table(x)))
  expect_null(reliability_profile(x)$estimate)
  expect_rp_error(reliability_draws(metadata, mean = rep(0, 4)), class = "rp_scope_error",
                  regexp = "variance-based")
  expect_null(.rp_brms_parts(brms_metadata(y | se(s) ~ 0 + x + (1 | site)))$mean_variable)
  expect_null(.rp_brms_parts(brms_metadata(y | se(s) ~ factor(x) + (1 | site)))$mean_variable)
  expect_null(.rp_brms_parts(brms_metadata(y | se(s) ~ stats::poly(x, 2) + (1 | site)))$mean_variable)
})

test_that("brms formula objects are read like plain formulas", {
  skip_if_not_installed("brms")
  x <- brms_metadata()
  x$formula <- brms::bf(y | se(s) ~ 1 + (1 | site))
  x$family <- brms::brmsfamily("gaussian")
  expect_identical(.rp_brms_parts(x)$sd_variable, "sd_site__Intercept")
  x$formula <- brms::bf(y | se(s) ~ 1 + x + (1 | site))
  expect_null(.rp_brms_parts(x)$mean_variable)
  x$formula <- brms::bf(y | se(s) ~ 1 + (1 | site), sigma ~ x)
  expect_rp_error(.rp_brms_parts(x), class = "rp_scope_error")
  x$formula <- brms::bf(y | se(s) ~ 1 + (1 | site) + ar(time = x))
  expect_rp_error(.rp_brms_parts(x), class = "rp_scope_error", regexp = "Autocorrelation")
})

test_that("a common mean written as 0 + Intercept is read; a mean fixed at 0 and prior draws are refused", {
  expect_identical(.rp_brms_parts(brms_metadata(y | se(s) ~ 0 + Intercept + (1 | site)))$mean_variable,
                   "b_Intercept")
  expect_rp_error(.rp_brms_parts(brms_metadata(y | se(s) ~ 0 + (1 | site))), "rp_scope_error",
                  regexp = "intercept")
  prior_only <- brms_metadata()
  prior_only$prior <- structure(data.frame(), sample_prior = "only")
  expect_rp_error(.rp_brms_parts(prior_only), "rp_scope_error", regexp = "prior")
})

test_that("a formula with an offset but no intercept is refused as having no intercept", {
  expect_rp_error(.rp_brms_parts(brms_metadata(y | se(s) ~ 0 + offset(x) + (1 | site))),
                  "rp_scope_error", regexp = "no intercept or coefficient")
})
