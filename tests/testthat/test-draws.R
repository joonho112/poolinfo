test_that("draws of sigma_tau give the diagnostics at each draw", {
  withr::local_seed(131)
  s <- abs(rnorm(3000, 8, 4))
  rd <- reliability_draws(s, se = eight_schools$se, model_label = "a test")
  expect_s3_class(rd, "reliability_draws")
  k <- .rp_draws_kernel(eight_schools$se^2)
  f <- .rp_functionals_at(k, s)
  ps <- profile_summary(rd)
  expect_equal(ps$mean[ps$quantity == "I"], mean(f[, "I"]))
  expect_equal(ps$mean[ps$quantity == "H"], mean(f[, "H"]))
  expect_equal(ps$lower[ps$quantity == "I"], unname(stats::quantile(f[, "I"], 0.025)))
  expect_true(all(ps$mcse[ps$quantity %in% c("I", "mean", "H", "J_eff")] > 0))
  expect_equal(prob_above(rd), mean(f[, "I"] > 0.5), ignore_attr = TRUE)
  expect_gt(attr(prob_above(rd), "mcse"), 0)
  expect_snapshot(print(rd, use_cli = FALSE))
})

test_that("draws from the exact posterior agree with its integrals", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  d <- profile_draws(fit, ndraws = 20000, seed = 7)
  rd <- reliability_draws(d$sigma_tau, se = eight_schools$se, estimate = eight_schools$estimate,
                          mean = d$effect_mean)
  for (q in c("I", "mean", "H", "J_eff")) {
    est <- rd$summary[rd$summary$quantity == q, ]
    exact <- fit$summary$mean[fit$summary$quantity == q]
    expect_lt(abs(est$mean - exact), 4.5 * est$mcse)
  }
  tab <- site_table(rd)
  ref <- fit$sites
  se_mc <- ref$effect_sd / sqrt(20000)
  expect_true(all(abs(tab$effect - ref$effect) < 5 * se_mc))
  expect_equal(tab$effect_sd, ref$effect_sd, tolerance = 0.03)
})

test_that("draws objects of the posterior package are read by variable name", {
  skip_if_not_installed("posterior")
  withr::local_seed(132)
  arr <- posterior::draws_array(sd_school = abs(stats::rnorm(800, 8, 3)), mu = stats::rnorm(800, 8, 5),
                                .nchains = 4)
  rd <- reliability_draws(arr, se = eight_schools$se, variable = "sd_school", mean = "mu",
                          estimate = eight_schools$estimate)
  expect_identical(nrow(rd$draws), 800L)
  expect_identical(length(unique(rd$chains)), 4L)
  # The one variable that looks like a standard deviation is found without a name.
  rd2 <- reliability_draws(arr, se = eight_schools$se)
  expect_equal(rd2$draws$sigma_tau, rd$draws$sigma_tau)
  expect_rp_error(reliability_draws(arr, se = eight_schools$se, variable = "tau_sd"),
                  regexp = "no variable `tau_sd`")
})

test_that("reliability_draws() checks its draws", {
  expect_rp_error(reliability_draws(c(1, -1), se = c(1, 1)), regexp = "cannot be negative")
  expect_rp_error(reliability_draws(c(1, 2), se = c(1, 1), mean = 1), regexp = "one draw for each")
  expect_rp_error(reliability_draws(c(1, NA), se = c(1, 1)), regexp = "finite")
  expect_rp_error(reliability_draws(5, se = c(1, 1)), regexp = "At least 2 draws")
  skip_if_not_installed("posterior")
  negative <- posterior::draws_df(tau = c(0.5, -0.2, 0.3))
  expect_rp_error(reliability_draws(negative, se = c(1, 1)), regexp = "cannot be negative")
})

test_that("draws with equal standard errors have a dispersion of 0", {
  d <- reliability_draws(c(0.5, 1, 2, 4), se = rep(0.3, 5))
  expect_identical(d$summary$mean[d$summary$quantity == "H"], 0)
})
