# Posteriors whose scale is far from that of the standard errors, units far
# from 1, and estimates far from 0.

test_that("a prior scale far below the standard errors gives the brute-force posterior", {
  withr::local_seed(201)
  v <- exp(rnorm(10, 0, 0.5))
  y <- rnorm(10, 0, sqrt(1 + v))
  for (s in c(1e-5, 1e-7)) {
    prior <- sd_prior("half_t", scale = s, df = 3)
    ref <- brute_rule(y, v, prior, lo = s * 1e-3, hi = 1e6, panels = 6000L)
    ref_sigma <- brute_mean(ref, identity)
    for (engine in c("nodes", "adaptive")) {
      fit <- reliability_posterior(y, se2 = v, prior = prior, control = posterior_control(engine = engine))
      expect_true(all(fit$checks$status == "ok"))
      expect_lt(abs(summary_value(fit, "sigma_tau") / ref_sigma - 1), 1e-7)
    }
  }
})

test_that("a uniform prior with its upper limit far below the standard errors gives a fit", {
  withr::local_seed(202)
  v <- exp(rnorm(10, 0, 0.5))
  y <- rnorm(10, 0, sqrt(1 + v))
  for (upper in c(3e-4, 1e-6)) {
    prior <- sd_prior("uniform", upper = upper)
    nodes <- reliability_posterior(y, se2 = v, prior = prior)
    adaptive <- reliability_posterior(y, se2 = v, prior = prior,
                                      control = posterior_control(engine = "adaptive"))
    expect_true(all(nodes$checks$status == "ok"))
    # The likelihood hardly changes over [0, upper]: the posterior is nearly
    # the uniform prior.
    expect_lt(abs(summary_value(nodes, "sigma_tau") / (upper / 2) - 1), 1e-5)
    expect_lt(abs(summary_value(adaptive, "sigma_tau") / summary_value(nodes, "sigma_tau") - 1), 1e-9)
  }
})

test_that("two sites far more precise than the others give the brute-force posterior", {
  withr::local_seed(203)
  yb <- rnorm(8, 0, 1.5)
  for (se_min in c(1e-4, 1e-6)) {
    v <- c(se_min^2, se_min^2, rep(1, 8))
    y <- c(0.3, 0.3 + 0.2 * se_min, 0.3 + yb)
    ref <- brute_rule(y, v, lo = se_min * 1e-3, hi = 1e6, panels = 8000L)
    I_at <- function(s) .rp_functionals_at(ref$kernel, s)[, "I"]
    fit <- reliability_posterior(y, se2 = v)
    expect_true(all(fit$checks$status == "ok"))
    expect_lt(abs(summary_value(fit, "I") - brute_mean(ref, I_at)), 1e-9)
    expect_lt(abs(brute_cdf(ref, fit$sigma_quantiles[[1L]]) - 0.025), 1e-9)
  }
})

test_that("estimates far from 0 keep their standard deviations and intervals", {
  y <- eight_schools$estimate
  se <- eight_schools$se
  base <- reliability_posterior(y, se = se)
  base_sites <- site_table(base)
  for (shift in c(1e8, 1e10)) {
    fit <- reliability_posterior(y + shift, se = se)
    expect_equal(summary_value(fit, "effect_mean", "sd"), summary_value(base, "effect_mean", "sd"),
                 tolerance = 1e-6)
    sites <- site_table(fit)
    expect_equal(sites$effect_sd, base_sites$effect_sd, tolerance = 1e-6)
    expect_equal(sites$effect_lower - shift, base_sites$effect_lower, tolerance = 1e-5)
    expect_equal(sites$effect_upper - shift, base_sites$effect_upper, tolerance = 1e-5)
  }
})

test_that("units from 1e-100 to 1e100 give the same diagnostics with both engines", {
  y <- eight_schools$estimate
  se <- eight_schools$se
  base <- reliability_posterior(y, se = se)
  rows <- c("I", "mean", "harmonic", "H", "J_eff")
  for (unit in c(1e-100, 1e-10, 1e10, 1e100)) {
    for (engine in c("nodes", "adaptive")) {
      fit <- expect_no_warning(reliability_posterior(y * unit, se = se * unit,
                                                     control = posterior_control(engine = engine)))
      expect_equal(fit$summary$mean[fit$summary$quantity %in% rows],
                   base$summary$mean[base$summary$quantity %in% rows], tolerance = 1e-9)
      expect_equal(summary_value(fit, "sigma_tau") / unit, summary_value(base, "sigma_tau"), tolerance = 1e-9)
      expect_equal(summary_value(fit, "effect_mean", "sd") / unit, summary_value(base, "effect_mean", "sd"),
                   tolerance = 1e-8)
      expect_equal(fit$sigma_quantiles / unit, base$sigma_quantiles, tolerance = 1e-9)
    }
  }
})

test_that("sampling variances over 300 orders of magnitude give a fit", {
  withr::local_seed(204)
  v <- c(1, 1e100, 1e200, 1e300)
  fit <- reliability_posterior(rnorm(4, 0, sqrt(v)), se2 = v)
  s <- fit$summary
  expect_true(all(is.finite(s$mean[s$quantity %in% c("I", "mean", "harmonic", "H", "J_eff")])))
  expect_true(is.finite(summary_value(fit, "effect_mean")))
})

test_that("a covariate on a large scale gives the interval of its coefficient", {
  withr::local_seed(205)
  x0 <- rnorm(10)
  v <- exp(rnorm(10, 0, 0.5))
  y <- 0.3 + 0.5 * x0 + rnorm(10, 0, sqrt(0.3 + v))
  base <- reliability_posterior(y, se2 = v, data = data.frame(x = x0), mean_model = ~x)
  for (unit in c(1e10, 1e13)) {
    fit <- reliability_posterior(y, se2 = v, data = data.frame(x = x0 * unit), mean_model = ~x)
    expect_equal(summary_value(fit, "coef_x", "lower") * unit, summary_value(base, "coef_x", "lower"),
                 tolerance = 1e-9)
    expect_equal(summary_value(fit, "coef_x", "upper") * unit, summary_value(base, "coef_x", "upper"),
                 tolerance = 1e-9)
  }
})

test_that("a prior with a heavy tail gives the same mean of sigma_tau with both engines", {
  prior <- sd_prior("half_t", scale = 1, df = 0.3)
  y <- c(0.3, -0.2)
  v <- c(0.04, 0.09)
  nodes <- suppressMessages(reliability_posterior(y, se2 = v, prior = prior))
  adaptive <- suppressMessages(reliability_posterior(y, se2 = v, prior = prior,
                                                     control = posterior_control(engine = "adaptive")))
  expect_equal(summary_value(adaptive, "sigma_tau"), summary_value(nodes, "sigma_tau"), tolerance = 1e-8)
  expect_true(all(adaptive$checks$status == "ok"))
})
