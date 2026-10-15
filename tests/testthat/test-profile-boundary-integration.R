test_that("external site tables retain subnormal reliabilities and effects", {
  sigma <- exp(-360)
  v <- c(1, 2)
  expected <- exp(2 * log(sigma) - log(v))
  x <- reliability_draws(rep(sigma, 4), se2 = v, estimate = c(1e100, 2e100), mean = rep(0, 4))
  sites <- site_table(x)
  expect_true(all(expected > 0))
  expect_identical(sites$reliability, expected)
  expect_identical(sites$reliability_lower, expected)
  expect_identical(sites$reliability_upper, expected)
  expect_equal(sites$effect / (expected * c(1e100, 2e100)), c(1, 1), tolerance = 1e-12)
  expect_true(all(sites$effect_sd > 0))
  expect_equal(sites$effect_sd^2 / (v * expected), c(1, 1), tolerance = 1e-8)
})

test_that("profiles at sigma_tau = 0 handle tiny sampling variances", {
  v <- c(1e-320, 2e-320)
  profile <- suppressMessages(reliability_profile(se2 = v, sigma_tau = 0, estimate = c(-1, 2)))
  sites <- site_table(profile)
  expect_equal(sites$share, c(2 / 3, 1 / 3), tolerance = 1e-12)
  expect_equal(sites$effect, c(0, 0), tolerance = 1e-12)
  expect_true(all(is.finite(as.matrix(sites[-1]))))
  expect_true(all(sites$effect_sd > 0))
  # The conditional mean variance is (1/v1 + 1/v2)^(-1), or 2*v1/3.
  # Compare after rescaling because subnormal variance has only a few digits.
  expect_equal(sites$effect_sd / sqrt(v[[1L]]), rep(sqrt(2 / 3), 2), tolerance = 1e-3)
})

test_that("site effects at a huge sigma_tau stay finite", {
  profile <- reliability_profile(se2 = c(1, 4), sigma_tau = 1e200, estimate = c(-1, 2))
  sites <- site_table(profile)
  expect_true(all(is.finite(as.matrix(sites[-1]))))
  expect_equal(sites$effect, c(-1, 2), tolerance = 1e-12)
  expect_equal(sites$effect_sd, c(1, 2), tolerance = 1e-12)
})
