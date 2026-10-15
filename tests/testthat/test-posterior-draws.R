test_that("draws are reproducible and leave the session's generator as it was", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  set.seed(5)
  before <- .Random.seed
  a <- profile_draws(fit, ndraws = 500, seed = 11)
  expect_identical(.Random.seed, before)
  b <- profile_draws(fit, ndraws = 500, seed = 11)
  expect_identical(a, b)
  expect_identical(attr(a, "seed"), 11L)
  expect_named(a, c("draw", "sigma_tau", "effect_mean", "I", "mean", "H", "J_eff"))
  c1 <- profile_draws(fit, ndraws = 10)
  expect_true(is.integer(attr(c1, "seed")))
})

test_that("draws agree with the integrals within their Monte Carlo error", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  d <- profile_draws(fit, ndraws = 20000, seed = 3, sites = TRUE)
  for (q in c("I", "mean", "H", "J_eff")) {
    exact <- fit$summary$mean[fit$summary$quantity == q]
    mcse <- stats::sd(d[[q]]) / sqrt(nrow(d))
    expect_lt(abs(mean(d[[q]]) - exact), 4.5 * mcse)
  }
  effects <- d[, paste0("effect_", 1:8)]
  mcse <- apply(effects, 2, stats::sd) / sqrt(nrow(d))
  expect_true(all(abs(colMeans(effects) - fit$sites$effect) < 4.5 * mcse))
  # The inverse of the table of the CDF is within 2e-8 of the exact CDF.
  expect_lt(attr(d, "inverse_cdf_error"), 2e-8)
  # Draws of sigma_tau have the right median.
  p_med <- mean(d$sigma_tau <= fit$sigma_quantiles[[2]])
  expect_lt(abs(p_med - 0.5), 4.5 * sqrt(0.25 / nrow(d)))
})

test_that("draws from a fit of the adaptive engine are those of the node rule", {
  nodes <- reliability_posterior(estimate, se = se, data = eight_schools)
  adaptive <- reliability_posterior(estimate, se = se, data = eight_schools,
                                    control = posterior_control(engine = "adaptive"))
  a <- profile_draws(adaptive, ndraws = 2000, seed = 4)
  b <- profile_draws(nodes, ndraws = 2000, seed = 4)
  expect_lt(max(abs(a$sigma_tau / b$sigma_tau - 1)), 1e-6)
  law <- .rp_fit_law(adaptive)
  expect_equal(law$quantile(c(1e-12, 1e-6)), .rp_fit_law(nodes)$quantile(c(1e-12, 1e-6)), tolerance = 1e-6)
})

test_that("site columns are named after the sites", {
  fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
  d <- profile_draws(fit, ndraws = 5, seed = 1, sites = TRUE)
  expect_true(all(paste0("effect_", LETTERS[1:8]) %in% names(d)))
  expect_rp_error(profile_draws(fit, ndraws = 0), regexp = "at least 1")
  expect_rp_error(profile_draws("x"), regexp = "reliability_posterior")
})
