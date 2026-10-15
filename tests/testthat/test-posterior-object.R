fit8 <- function(...) {
  reliability_posterior(estimate, se = se, site = school, data = eight_schools, ...)
}

test_that("a posterior fit has its parts and prints its words", {
  fit <- fit8()
  expect_s3_class(fit, "reliability_posterior")
  expect_identical(fit$summary$quantity,
                   c("I", "mean", "harmonic", "H", "J_eff", "V", "sigma_tau", "sigma_tau2", "effect_mean"))
  expect_identical(fit$sites$site, LETTERS[1:8])
  expect_true(all(fit$checks$status == "ok"))
  expect_snapshot(print(fit, use_cli = FALSE))
  expect_snapshot(print(summary(fit), use_cli = FALSE))
})

test_that("the data can be given as vectors or as columns", {
  a <- fit8()
  b <- reliability_posterior(eight_schools$estimate, se = eight_schools$se, site = eight_schools$school)
  expect_equal(a$summary, b$summary)
  c2 <- reliability_posterior(estimate, se2 = se^2, data = eight_schools)
  expect_equal(a$summary$mean, c2$summary$mean, tolerance = 1e-12)
  expect_identical(c2$sites$site, as.character(1:8))
  expect_rp_error(reliability_posterior(estimate, se = sd_col, data = eight_schools),
                  regexp = "could not be found")
  expect_rp_error(reliability_posterior(se = eight_schools$se), regexp = "`estimate` is missing")
  expect_rp_error(reliability_posterior(1:3, se = 1:2), regexp = "one value for each of the 3")
  expect_rp_error(reliability_posterior(1:3, se = 1:3, data = list(a = 1)), regexp = "data frame")
  expect_rp_error(fit8(level = 1), regexp = "between 0 and 1")
  expect_rp_error(fit8(prior = "flat"), regexp = "sd_prior")
  expect_rp_error(fit8(control = list(engine = "nodes")), regexp = "posterior_control")
})

test_that("missing values fail by default and are left out on request", {
  d <- eight_schools
  d$estimate[3] <- NA
  expect_rp_error(reliability_posterior(estimate, se = se, data = d), regexp = "missing value")
  expect_message(fit <- reliability_posterior(estimate, se = se, site = school, data = d, na_action = "omit"),
                 class = "rp_message_omitted")
  expect_identical(fit$sites$site, LETTERS[1:8][-3])
  expect_identical(fit$data$omitted, 1L)
})

test_that("correlated estimates are outside the scope", {
  expect_rp_error(reliability_posterior(1:3, se2 = diag(3)), "rp_scope_error")
})

test_that("the flat prior needs 3 sites, and the moments of sigma_tau need more", {
  err <- expect_rp_error(reliability_posterior(c(1, 2), se = c(1, 1)), regexp = "at least 3 sites")
  expect_match(conditionMessage(err), "sd_prior", fixed = TRUE)
  fit3 <- reliability_posterior(c(1, 5, -2), se = c(1, 2, 1.5))
  s <- fit3$summary
  expect_true(is.na(s$mean[s$quantity == "sigma_tau"]))
  expect_true(is.na(s$mean[s$quantity == "sigma_tau2"]))
  expect_false(anyNA(s$mean[s$quantity %in% c("I", "mean", "H", "J_eff")]))
  expect_false(anyNA(s$lower))
  expect_true(is.na(fit3$plugin[["I"]]))
  expect_message(pg <- plugin_gap(fit3), class = "rp_message_no_moment")
  expect_true(all(is.na(pg$at_mean_sigma2)))
  expect_rp_error(reliability_profile(fit3), regexp = "does not exist")
  fit4 <- reliability_posterior(c(1, 5, -2, 3), se = c(1, 2, 1.5, 1))
  expect_false(is.na(fit4$summary$mean[fit4$summary$quantity == "sigma_tau"]))
  expect_true(is.na(fit4$summary$mean[fit4$summary$quantity == "sigma_tau2"]))
  fit5 <- reliability_posterior(c(1, 5, -2, 3, 0), se = c(1, 2, 1.5, 1, 1))
  expect_false(is.na(fit5$summary$mean[fit5$summary$quantity == "sigma_tau2"]))
  expect_true(is.na(fit5$summary$sd[fit5$summary$quantity == "sigma_tau2"]))
  expect_identical(fit5$moments$exists, c(TRUE, TRUE, FALSE))
})

test_that("the posterior mean of I is below I at the posterior mean of sigma_tau^2 (Proposition 3)", {
  fit <- fit8()
  pg <- plugin_gap(fit)
  expect_gt(pg$difference[pg$functional == "I"], 0)
  expect_equal(attr(pg, "sigma2"), fit$summary$mean[fit$summary$quantity == "sigma_tau2"])
})

test_that("probabilities agree with the quantiles of sigma_tau", {
  fit <- fit8()
  q <- fit$sigma_quantiles
  # The probability that sigma_tau exceeds its 97.5% quantile.
  expect_equal(prob_above(fit, q[[3L]], "sigma_tau"), 0.025, tolerance = 1e-9)
  I_975 <- fit$quantiles$upper[fit$quantiles$quantity == "I"]
  expect_equal(prob_above(fit, I_975), 0.025, tolerance = 1e-9)
  mean_975 <- fit$quantiles$upper[fit$quantiles$quantity == "mean"]
  expect_equal(prob_above(fit, mean_975, "mean"), 0.025, tolerance = 1e-8)
  jeff_med <- fit$quantiles$median[fit$quantiles$quantity == "J_eff"]
  expect_equal(prob_above(fit, jeff_med, "J_eff"), 0.5, tolerance = 1e-8)
  expect_identical(prob_above(fit), fit$probability[["I_above"]])
  expect_rp_error(prob_above(fit, which = "mean"), regexp = "`threshold` is needed")
  expect_rp_error(prob_above(fit, 20, "J_eff"), regexp = "between 1 and 8")
  expect_rp_error(prob_above(1), regexp = "no method")
})

test_that("summaries at another level are recomputed", {
  fit <- fit8()
  s90 <- profile_summary(fit, level = 0.9)
  expect_identical(attr(s90, "level"), 0.9)
  expect_equal(s90$mean, fit$summary$mean)
  refit <- fit8(level = 0.9)
  expect_equal(s90$lower, refit$summary$lower, tolerance = 1e-10)
  expect_equal(s90$upper, refit$summary$upper, tolerance = 1e-10)
  expect_equal(profile_summary(fit)[, 1:5], fit$summary, ignore_attr = TRUE)
})

test_that("update() refits with another prior or level", {
  fit <- fit8()
  up <- update(fit, prior = sd_prior("half_normal", scale = 10))
  expect_identical(up$prior$type, "half_normal")
  expect_false(isTRUE(all.equal(up$summary$mean, fit$summary$mean)))
  expect_identical(update(fit, level = 0.8)$level, 0.8)
})

test_that("the two engines agree on eight schools", {
  a <- fit8(control = posterior_control(engine = "adaptive"))
  n <- fit8()
  rows <- c("I", "mean", "harmonic", "H")
  expect_equal(a$summary$mean[a$summary$quantity %in% rows], n$summary$mean[n$summary$quantity %in% rows],
               tolerance = tol_bounded)
  expect_lt(max(abs(a$sites$effect - n$sites$effect)), tol_effect_scale(eight_schools$se^2))
  expect_lt(max(abs(a$sites$kl_term - n$sites$kl_term)), tol_bounded)
})

test_that("check_numerics() returns the table of checks", {
  ch <- check_numerics(fit8())
  expect_named(ch, c("check", "value", "limit", "status"))
  expect_true(all(ch$status == "ok"))
  expect_rp_error(check_numerics(1), regexp = "reliability_posterior")
})
