test_that("the densities of the diagnostics give the mass between their quantiles", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  q <- fit$quantiles
  for (which in c("I", "mean", "J_eff", "sigma_tau", "H")) {
    d <- functional_density(fit, which, n = 4001)
    expect_named(d, c("value", "density"))
    expect_true(all(d$density >= 0))
    # Between the median and the 97.5% quantile the density is smooth; near
    # the lower end it can grow without bound, where the trapezoid rule would
    # overstate the mass.
    lo <- q$median[q$quantity == which]
    hi <- q$upper[q$quantity == which]
    grid <- seq(lo, hi, length.out = 2001)
    dens <- stats::approx(d$value, d$density, xout = grid)$y
    mass <- sum(diff(grid) * (dens[-1] + dens[-2001]) / 2)
    expect_equal(mass, 0.475, tolerance = 1e-4, label = which)
  }
})

test_that("the density of I is the derivative of its distribution function", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  d <- functional_density(fit, "I", n = 50)
  k <- .rp_fit_kernel(fit)
  law <- .rp_fit_law(fit, k)
  cdf_I <- function(i) law$cdf(sqrt(exp(k$log_gm) * i / (1 - i)))
  h <- 1e-6
  at <- d$value[c(10, 25, 40)]
  numeric_density <- (cdf_I(at + h) - cdf_I(at - h)) / (2 * h)
  expect_equal(d$density[c(10, 25, 40)], numeric_density, tolerance = 1e-5)
})

test_that("the density of H adds its two branches", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  d <- functional_density(fit, "H", n = 40)
  k <- .rp_fit_kernel(fit)
  law <- .rp_fit_law(fit, k)
  br <- .rp_H_branches(k, .rp_turning_points(k$d))
  h <- 1e-7
  at <- d$value[c(8, 20, 32)]
  numeric_density <- (.rp_H_cdf(law, k, br, at + h) - .rp_H_cdf(law, k, br, at - h)) / (2 * h)
  expect_equal(d$density[c(8, 20, 32)], numeric_density, tolerance = 1e-4)
})

test_that("functional_density() checks its arguments", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  expect_rp_error(functional_density(fit, "V"), regexp = "one of")
  expect_rp_error(functional_density(fit, n = 2), regexp = "at least 10")
  expect_rp_error(functional_density(1), regexp = "reliability_posterior")
  equal <- reliability_posterior(c(1, 2, 3, 4), se = rep(1, 4))
  expect_rp_error(functional_density(equal, "H"), regexp = "no density")
})
