profile_model_data <- function() {
  x <- seq(-2, 2, length.out = 12)
  data.frame(y = 2 * x + c(-.05, .12, -.07, .03, -.08, .04, .06, -.03, .09, -.10, .01, .02),
             v = seq(.02, .06, length.out = 12), x = x)
}

fixed_model_reference <- function(y, v, X, sigma_tau) {
  w <- 1 / (v + sigma_tau^2)
  cov <- solve(crossprod(X, w * X))
  predicted <- as.vector(X %*% cov %*% crossprod(X, w * y))
  prediction_var <- unname(rowSums((X %*% cov) * X))
  S <- sigma_tau^2 / (sigma_tau^2 + v)
  list(mean = S * y + (1 - S) * predicted,
       var = v * S + (1 - S)^2 * prediction_var)
}

test_that("a single metafor moderator without an intercept keeps its mean model", {
  skip_if_not_installed("metafor")
  d <- profile_model_data()
  fit <- metafor::rma(y, v, mods = ~ 0 + x, data = d)
  actual <- reliability_posterior(fit)
  direct <- reliability_posterior(y, se2 = v, mean_model = ~ 0 + x, data = d)
  expect_equal(actual$data$X, fit$X, ignore_attr = TRUE)
  expect_equal(actual$summary, direct$summary, tolerance = 1e-10)
  expect_equal(site_table(actual), site_table(direct), tolerance = 1e-10)
  expect_lt(actual$summary$mean[actual$summary$quantity == "I"], .3)
  common <- metafor::rma(y, v, data = d)
  expect_null(reliability_posterior(common)$data$X)
  expect_equal(reliability_posterior(common)$summary,
               reliability_posterior(y, se2 = v, data = d)$summary)
})

test_that("fixed profiles of metafor mean models retain conditional site moments", {
  skip_if_not_installed("metafor")
  d <- profile_model_data()
  for (formula in list(~ x, ~ 0 + x)) {
    for (tau2 in c(0, .09)) {
      fit <- metafor::rma(y, v, mods = formula, data = d, tau2 = tau2)
      profile <- suppressMessages(reliability_profile(fit))
      reference <- fixed_model_reference(d$y, d$v, fit$X, sqrt(tau2))
      sites <- site_table(profile)
      expect_equal(profile$X, fit$X)
      expect_equal(sites$effect, reference$mean, tolerance = 1e-12)
      expect_equal(sites$effect_sd^2, reference$var, tolerance = 1e-12)
      plain <- suppressMessages(reliability_profile(se2 = d$v, sigma_tau = sqrt(tau2)))
      expect_identical(profile$stats, plain$stats)
      updated <- update(profile, sigma_tau = .4)
      updated_reference <- fixed_model_reference(d$y, d$v, fit$X, .4)
      expect_identical(updated$X, profile$X)
      expect_equal(site_table(updated)$effect, updated_reference$mean, tolerance = 1e-12)
      expect_equal(site_table(updated)$effect_sd^2, updated_reference$var, tolerance = 1e-12)
    }
  }
})

test_that("a profile from a posterior keeps the mean model in its site effects and pooling plot", {
  d <- profile_model_data()
  for (formula in list(~ x, ~ 0 + x)) {
    fit <- reliability_posterior(y, se2 = v, mean_model = formula, data = d)
    profile <- reliability_profile(fit)
    reference <- fixed_model_reference(d$y, d$v, fit$data$X, profile$sigma_tau)
    expect_identical(profile$X, fit$data$X)
    expect_equal(site_table(profile)$effect, reference$mean, tolerance = 1e-12)
    expect_equal(site_table(profile)$effect_sd^2, reference$var, tolerance = 1e-12)
    expect_equal(profile$stats,
                 reliability_profile(se2 = d$v, sigma_tau = profile$sigma_tau)$stats)
    plot <- plot_pooling(profile, order = "input", level = .8)
    expect_equal(plot$data$effect, rev(reference$mean), tolerance = 1e-12)
    expect_equal(plot$data$effect_lower,
                 rev(reference$mean - qnorm(.9) * sqrt(reference$var)), tolerance = 1e-12)
    expect_false(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomVline"), FALSE)))
  }
})

test_that("fixed and equal effect metafor fits have an explicit scope error", {
  skip_if_not_installed("metafor")
  d <- profile_model_data()
  for (method in c("FE", "EE", "CE")) {
    fit <- metafor::rma(y, v, data = d, method = method)
    expect_rp_error(reliability_profile(fit), "rp_scope_error", regexp = "random-effects")
    expect_rp_error(reliability_posterior(fit), "rp_scope_error", regexp = "random-effects")
    expect_rp_error(map_tau2_interval(fit), "rp_scope_error", regexp = "random-effects")
    expect_s3_class(reliability_posterior(fit$yi, se2 = fit$vi), "reliability_posterior")
  }
})
