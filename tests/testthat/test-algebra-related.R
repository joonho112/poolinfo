test_that("the typical variance is (J - 1) W / (W^2 - W_2)", {
  withr::local_seed(51)
  for (i in 1:100) {
    v <- exp(runif(sample(2:50, 1), -6, 2))
    w <- 1 / v
    expect_equal(typical_variance(se2 = v), (length(v) - 1) * sum(w) / (sum(w)^2 - sum(w^2)),
                 tolerance = 1e-10)
    expect_gte(typical_variance(se2 = v), 1 / mean(1 / v) * (1 - 1e-12))
  }
  expect_equal(typical_variance(se2 = rep(0.04, 6)), 0.04, tolerance = 1e-15)
  expect_equal(typical_variance(se = c(0.1, 0.2)), typical_variance(se2 = c(0.01, 0.04)))
})

test_that("the typical variance keeps its digits when one weight dominates", {
  v <- c(1e-10, rep(1, 9))
  w <- 1 / v
  exact <- 9 * sum(w) / (2 * w[[1]] * 9 + 2 * choose(9, 2))
  expect_equal(typical_variance(se2 = v), exact, tolerance = 1e-14)
})

test_that("typical variance survives extreme common scales and underflowed weight ratios", {
  # For two sites v_typ is their arithmetic mean, independent of inverse
  # variance weights. Ratios make zero/Inf failures visible at every scale.
  cases <- list(c(1e308, 1e308), c(1e-100, 1e200),
                c(1e-300, 1e300), c(1e-320, 1))
  for (v in cases) {
    actual <- typical_variance(se2 = v)
    expected <- mean(v)
    expect_true(is.finite(actual) && actual > 0)
    expect_equal(actual / expected, 1, tolerance = 1e-12)
  }
  expect_identical(typical_variance(se2 = rep(1e308, 2)), 1e308)
  # For one tiny variance and two equal large variances the limit is v/2.
  expect_equal(typical_variance(se2 = c(1e-300, 1e300, 1e300)) / 5e299,
               1, tolerance = 1e-12)
})

test_that("metafor's I^2 is the reliability at the typical variance", {
  skip_if_not_installed("metafor")
  withr::local_seed(52)
  vi <- runif(15, 0.01, 0.2)
  yi <- rnorm(15, 0.2, sqrt(vi + 0.05))
  fit <- metafor::rma(yi = yi, vi = vi, method = "REML")
  vt <- typical_variance(se2 = vi)
  expect_equal(fit$I2 / 100, fit$tau2 / (fit$tau2 + vt), tolerance = 1e-10)
})

test_that("the typical variance needs two sites", {
  expect_rp_error(typical_variance(se = 0.1), regexp = "at least 2 numbers")
})

test_that("the benchmark risk is sigma_tau^2 times one minus the mean reliability (Proposition 2(iii))", {
  withr::local_seed(53)
  for (i in 1:100) {
    p <- random_profile(log_ratio = 3, sigma_tau = exp(runif(1, -2, 1)))
    expect_equal(benchmark_risk(se2 = p$v, sigma_tau = p$sigma_tau),
                 p$sigma_tau^2 * (1 - mean(p$S)), tolerance = 1e-12)
  }
  expect_identical(benchmark_risk(se = c(0.1, 0.2), sigma_tau = 0), 0)
  expect_equal(benchmark_risk(se2 = c(0.01, 0.04), sigma_tau = 1e200), 0.025)
  expect_rp_error(benchmark_risk(se = 0.1), regexp = "`sigma_tau` is missing")
})

test_that("benchmark risks and reliabilities at the typical variance avoid intermediate overflow", {
  risk <- benchmark_risk(se2 = 1e200, sigma_tau = 1e-100)
  expect_gt(risk, 0)
  expect_equal(risk / 1e-200, 1, tolerance = 1e-12)
  p <- reliability_profile(se2 = c(1, 4), sigma_tau = 1e200)
  expect_identical(unname(p$stats[["R_typ"]]), 1)
  expect_equal(unname(p$stats[["benchmark_risk"]]) / 2.5, 1, tolerance = 1e-12)
  p_small <- reliability_profile(se2 = c(1e200, 1e200), sigma_tau = 1e-100)
  expect_equal(unname(p_small$stats[["benchmark_risk"]]) / 1e-200, 1, tolerance = 1e-12)
  # The finite ratio to v still matters on both sides of square overflow.
  expect_equal(benchmark_risk(se2 = 1e308, sigma_tau = 1e153) / (1e308 / 101),
               1, tolerance = 1e-12)
  expect_equal(benchmark_risk(se2 = 1e308, sigma_tau = 1e155) / (1e308 / 1.01),
               1, tolerance = 1e-12)
  p_zero <- suppressMessages(reliability_profile(se2 = c(1e-320, 1), sigma_tau = 0))
  expect_true(all(is.finite(p_zero$stats)))
  expect_identical(unname(p_zero$stats[["J_eff"]]), 1)
})

test_that("the Burg and squared dispersions have their closed forms", {
  withr::local_seed(54)
  S <- random_profile(J = 15)$S
  gm <- exp(mean(log(S)))
  hm <- 1 / mean(1 / S)
  expect_equal(bregman_dispersion(S, "burg", "first"), log(gm / hm), tolerance = 1e-12)
  expect_equal(bregman_dispersion(S, "burg", "second"), log(mean(S) / gm), tolerance = 1e-12)
  expect_equal(bregman_dispersion(S, "squared", "first"), mean((S - mean(S))^2) / 2)
  expect_equal(bregman_dispersion(S, "squared", "second"), mean((S - mean(S))^2) / 2)
})
