test_that("reliabilities are sigma_tau^2 / (sigma_tau^2 + v)", {
  se <- c(a = 0.1, b = 0.15, c = 0.2, d = 0.3)
  S <- site_reliability(se = se, sigma_tau = 0.2)
  expect_equal(S, 0.04 / (0.04 + se^2), tolerance = 1e-15)
  expect_named(S, names(se))
  expect_identical(site_reliability(se2 = se^2, sigma_tau = 0.2), S)
  expect_equal(unname(S), c(0.8, 0.64, 0.5, 0.04 / 0.13), tolerance = 1e-15)
})

test_that("a between-site SD of 0 gives zeros and a message", {
  expect_message(S <- site_reliability(se = c(0.1, 0.2), sigma_tau = 0),
                 class = "rp_message_sigma_zero")
  expect_identical(S, c(0, 0))
})

test_that("site_reliability() checks its arguments", {
  expect_rp_error(site_reliability(se = 0.1), regexp = "`sigma_tau` is missing")
  expect_rp_error(site_reliability(se = 0.1, sigma_tau = -1), regexp = "at least 0")
  expect_rp_error(site_reliability(se = 0.1, se2 = 0.01, sigma_tau = 1), regexp = "both supplied")
  expect_rp_error(site_reliability(sigma_tau = 1), regexp = "standard errors are missing")
  err <- expect_rp_error(site_reliability(se = c(0.1, 0), sigma_tau = 1))
  expect_identical(error_call_name(err), "site_reliability")
})

test_that("the logits are log(sigma^2) - log(v) and keep their digits near 0 and 1", {
  v <- c(1e-20, 1, 1e20)
  eta <- .rp_logit_from_variances(1, v)
  expect_equal(eta, -log(v), tolerance = 1e-15)
  # The reliability of the least precise site is 1e-20, which a logit taken
  # from S itself would also give, but 1 - S of the most precise site is
  # 1e-20, which S itself cannot hold.
  expect_equal(stats::plogis(-eta[[1L]]), 1e-20, tolerance = 1e-12)
})

test_that("centering gives deviations that sum to zero", {
  withr::local_seed(11)
  eta <- rnorm(50, 300, 5)
  p <- .rp_center(eta)
  expect_equal(p$m + p$d, eta, tolerance = 1e-15)
  expect_lt(abs(sum(p$d)), 1e-12)
  cv <- .rp_center_v(exp(eta / 100))
  expect_equal(cv$log_gm, mean(eta / 100), tolerance = 1e-15)
  expect_lt(abs(sum(cv$d)), 1e-12)
})

test_that("softplus is log(1 + exp(x)) without overflow", {
  x <- c(-800, -30, -1, 0, 1, 30, 800)
  expect_equal(.rp_softplus(x[2:6]), log1p(exp(x[2:6])), tolerance = 1e-15)
  expect_identical(.rp_softplus(800), 800)
  expect_equal(.rp_softplus(-800), exp(-800))
})

test_that("reliabilities keep values below the smallest normal number", {
  # exp(-736.827...) is representable; base plogis() returns zero here.
  S <- site_reliability(se2 = 1, sigma_tau = 1e-160)
  p <- reliability_profile(se2 = 1, sigma_tau = 1e-160)
  expected <- exp(2 * log(1e-160))
  expect_gt(expected, 0)
  expect_identical(S, expected)
  expect_identical(p$S, expected)
  expect_identical(unname(p$stats[["I"]]), expected)
  expect_identical(unname(p$stats[["trace"]]), expected)
})
