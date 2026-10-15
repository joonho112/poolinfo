# A few inputs at the edges, run with every test; test-engines-slow.R runs
# 400 of them, with up to 10,000 sites, when POOLINFO_SLOW is "true".

fit_is_sound <- function(fit) {
  s <- fit$summary
  bounded <- s$quantity %in% c("I", "mean", "harmonic", "H", "J_eff", "V")
  expect_true(all(is.finite(s$mean[bounded])))
  expect_true(all(is.finite(s$lower)) && all(is.finite(s$upper)))
  expect_true(all(s$lower <= s$upper))
  expect_gte(fit$probability[["I_above"]], 0)
  expect_lte(fit$probability[["I_above"]], 1)
  expect_true(all(is.finite(as.matrix(fit$sites[, -1]))))
  expect_true(all(fit$checks$status == "ok"))
}

test_that("variance ratios up to 1e8 and estimates far apart give sound fits", {
  withr::local_seed(101)
  for (ratio in c(1e2, 1e5, 1e8)) {
    v <- 0.01 * exp(seq(0, log(ratio), length.out = 12))
    y <- rnorm(12, 0, sqrt(0.02 + v))
    y[1] <- y[1] + 40 * sqrt(max(v))
    fit_is_sound(reliability_posterior(y, se2 = v))
  }
})

test_that("nearly equal variances give a dispersion of nearly 0", {
  withr::local_seed(102)
  for (spread in c(1e-12, 1e-8, 1e-4)) {
    v <- 0.04 * (1 + runif(10, -1, 1) * spread)
    fit <- reliability_posterior(rnorm(10, 0, 0.3), se2 = v)
    fit_is_sound(fit)
    expect_lt(fit$summary$upper[fit$summary$quantity == "H"], spread^2)
  }
  equal <- reliability_posterior(c(0.1, 0.4, -0.2, 0.3), se = rep(0.2, 4))
  expect_true(equal$h_branches$homogeneous)
  expect_identical(equal$summary$mean[equal$summary$quantity == "H"], 0)
})

test_that("estimates on very small and very large scales give the same diagnostics", {
  withr::local_seed(103)
  v <- runif(15, 0.5, 2)
  y <- rnorm(15, 0, sqrt(1 + v))
  base <- reliability_posterior(y, se2 = v)
  for (scale in c(1e-6, 1e6)) {
    fit <- reliability_posterior(y * scale, se2 = v * scale^2)
    fit_is_sound(fit)
    rows <- c("I", "mean", "harmonic", "H", "J_eff")
    expect_equal(fit$summary$mean[fit$summary$quantity %in% rows],
                 base$summary$mean[base$summary$quantity %in% rows], tolerance = 1e-9)
    expect_equal(fit$summary$mean[fit$summary$quantity == "sigma_tau"] / scale,
                 base$summary$mean[base$summary$quantity == "sigma_tau"], tolerance = 1e-9)
  }
})

test_that("data with no spread beyond the standard errors give a posterior near 0", {
  withr::local_seed(104)
  v <- runif(30, 0.5, 2)
  fit <- reliability_posterior(rep(0, 30), se2 = v)
  fit_is_sound(fit)
  expect_lt(fit$summary$mean[fit$summary$quantity == "I"], 0.2)
})
